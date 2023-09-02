clear all
clc
close all

% Code to compute pairwise CCGs for the new alert data

%% Load data structure here

datasetpaths{1} = 'E:\Data\H07\07-04-2016\PFC\Spontaneous\SUSpikesRS.mat';
datasetpaths{2} = 'E:\Data\A11\06-04-2016\PFC\Spontaneous\SUSpikesRS.mat';
tag = [{'spon'},{'anethspon'}];
[SU, n, dataInfo, utahMaps] = organize_MEAephysData(datasetpaths,'spon');
directories{1}.PFC = 'E:\Data\H07\07-04-2016\PFC';
directories{2}.PFC = 'E:\Data\A11\06-04-2016\PFC';

%% Load data structure here

% datasetpaths{1} = 'Z:\D98T\2010-06-01_09-21-03\spontaneous02\SUSpikesRS.mat';
% datasetpaths{2} = 'Y:\Makkay\2010-03-23_24\2010-03-24_22-05-27\spontaneous6\SUSpikesRS.mat';
% tag = [{'spon'},{'anethspon'}];
% [SU, n, dataInfo, utahMaps] = organize_MEAephysData(datasetpaths,'anethSpon');
% directories{1}.PFC = 'Z:\D98T\2010-06-01_09-21-03';
% directories{2}.PFC = 'Y:\Makkay\2010-03-23_24\2010-03-24_22-05-27';
%% Params

params.contamthreshold = 10000000000;
params.binsize = 0.025;
params.maxlag = 0.5;
%params.fs = 32556;
params.fs = 3e4;
params.datasets = 2;
rec = tag{1};
%% Start processing the CCGs

for i = 1:params.datasets
    
    nDataSet = i;
    msg=sprintf('Data set %d of %d starts...\n',i,nDataSet);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances and contamination values
    
    SU_PWdistances = SU(1,i).PWdistances; % distances
    
    SUcontam = utahMaps(1,i).SU.integratedTabel(:,5); % contamination indices
    
    SU_PWdistances = triu(SU_PWdistances); % Remove doubles
    
    % Identify pairs in distance bins
    
    dist = identifypairs(SU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    SUind = find(SUcontam<=params.contamthreshold);
    
    nunits = size(SU_PWdistances,1);
    
    % Eliminate low spiking cells
    
sc = combvec(1:nunits,1:nunits);
sc = sc';
csorted = sort(sc,2);
sitePairs = unique(csorted,'rows');

SUspikesByTimeT = SU(i).spikesByTime;
cell_pairs{i} = findValPairs3(dist,SU_PWdistances,SUind);

if strcmp(rec,'anethspon')

cd(directories{i}.PFC);
FieldSelectionEV = [1 1 1 0 1];
[TimeStamps, ~, Nttls, EventStrings] = Nlx2MatEV('Events.nev', FieldSelectionEV, 0, 1);
PFCevents.times = TimeStamps;
PFCevents.events = EventStrings;

%tagsNeeded = [{'sw21_movie02'}, {'sw21_movie02 end'}]; % INPUT THE REQUIRED TAG FROM THE ABOVE LIST HERE!! CHECK THE EXCEL LOGFILE BEFORE YOU INPUT THE TAGS
%tagsNeeded = [{'Movie1'}, {'Movie1 end'}];
if i == 1
    
tagsNeeded = [{'spontaneous02'}, {'spontaneous02 end'}]; % INPUT THE REQUIRED TAG FROM THE ABOVE LIST HERE!! CHECK THE EXCEL LOGFILE BEFORE YOU INPUT THE TAGS
else
    tagsNeeded = [{'spontaneous6'}, {'spontaneous6 end'}];
end
%tagsNeeded = [{'OriDisk1'}, {'OriDisk1 end'}];

%Use ismember to get indices of the epochs that did happen in the experiment
expEpochTags = cellstr(PFCevents.events);

[ispresent,epochLogical] = ismember(tagsNeeded,expEpochTags);

%epochLogical(epochLogical==0) = [];

epochinds = PFCevents.times(epochLogical); % THIS IS IN MICROSECONDS!!!

epochinds = epochinds/1e3; % in ms

corrsTrial{i} = sccgparnewanethspon(params,SUspikesByTimeT,sitePairs,epochinds);
else
%corrsTrial{i} = sccgparnewspon(params,SUspikesByTimeT,sitePairs,directories{i});
corrsTrial{i} = sccgparnewsponchopped(params,SUspikesByTimeT,sitePairs,directories{i});
end
end
%% Save Awake
cd('L:\projects\AbhilashD\Correlations\NewAwakeData')
save('AwakeSpontaneousCorrsChopped5Hz.mat','cell_pairs','corrsTrial','datasetpaths','dataInfo','utahMaps');

%% Save Anaesthetised

cd('L:\projects\AbhilashD\Correlations\NewAnethData')
save('AnethSpontaneousCorrs.mat','cell_pairs','corrsTrial','datasetpaths','dataInfo','utahMaps');

%% Plot Anaesthetised

dBins = 0.5:0.5:4;
lags = linspace(-0.5,0.5,501);

for i = 1:params.datasets
for j = 1:length(cell_pairs{i})
for k = 1:size(cell_pairs{i}{j},1)
pc = corrsTrial{i}.meanCCGs{cell_pairs{i}{j}(k,1),cell_pairs{i}{j}(k,2)};
pairCorrs(k,:) = pc;
end
corrsdBin(j,:) = nanmean(pairCorrs,1);
clear pairCorrs
end
dataSetBinCorrs{i} = corrsdBin;
end

figure(3)
subplot(2,2,[1 3])
for i = 1:8
hold on
plot(lags,normalise(dataSetBinCorrs{2}(i,:))+(i*1))
end
grid on
axis tight
box off
xlabel('lags in seconds')
ylabel('Proportionally increased coincidence/spk')
title('Makkay Task Synchrony by Distance')

subplot(2,2,[2 4])
for i = 1:8
    peaks(i) = max(dataSetBinCorrs{2}(i,:));
end
plot(dBins,peaks)
grid on
axis tight
box off
xlabel('distance in mm')
ylabel('CCG Peak (coincidences/spk')
title('Makkay Task Structure of Synchrony')

%% Plot Awake

dBins = 0.5:0.5:4;
lags = linspace(-0.5,0.5,501);

for i = 1:params.datasets
    
    for j = 1:length(cell_pairs{i}) % 8 distance bins
        
        for k = 1:size(cell_pairs{i}{j},1) % n Pairs in each distance bin
            
            if i == 1 && cell_pairs{i}{j}(k,1)==64 && cell_pairs{i}{j}(k,2)==65 || i == 1 && cell_pairs{i}{j}(k,1)==63 && cell_pairs{i}{j}(k,2)==66 || i == 2 && cell_pairs{i}{j}(k,1)==65 && cell_pairs{i}{j}(k,2)==67
            pairCorrs(k,:) = nan(1,501);
            else
            pairCorrs(k,:) = corrsTrial{i}.meanCCGs{cell_pairs{i}{j}(k,1),cell_pairs{i}{j}(k,2)}; % Collect the ccgs
            end
            
        end
        
        corrsdBin(j,:) = nanmean(pairCorrs,1);
        clear pairCorrs
        
    end
    
    dataSetBinCorrs{i} = corrsdBin;
    
end

fignames{1} = ['H07_07042016'];
fignames{2} = ['H07_06042016'];

for i = 1:2
figure(i)
subplot(2,2,[1 3])
for j = 1:8
hold on
plot(lags,normalise(dataSetBinCorrs{i}(j,:))+(j*1))
end
grid on
axis tight
box off
xlabel('lags in seconds')
ylabel('Proportionally increased coincidence/spk')
title([fignames{i} 'Spon Synchrony by Distance'])

subplot(2,2,[2 4])
for j = 1:8
    peaks(j) = trapz(dataSetBinCorrs{i}(j,:));
end
plot(dBins,peaks)
grid on
axis tight
box off
xlabel('distance in mm')
ylabel('CCG Peak (coincidences/spk')
title([fignames{i} 'Spon Structure of Synchrony'])

print('-dpng', '-r600', strcat([fignames{i} '1sChopped_5Hz_Spon_StructSynchrony_Integral'],'.png'));
end