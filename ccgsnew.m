clear all
clc
close all

% Code to compute pairwise CCGs for the new alert data

%% Load data structure here

datasetpaths{1} = 'E:\Data\H07\20161019\PFC\OriDisk1\SUSpikesByTime.mat';
datasetpaths{2} = 'E:\Data\H07\20161025\PFC\OriDisk1\SUSpikesByTime.mat';
datasetpaths{3} = 'E:\Data\A11\20161130\PFC\OriDisk1\SUSpikesByTime.mat';
datasetpaths{4} = 'E:\Data\A11\20161203\PFC\OriDisk1\SUSpikesByTime.mat';
[SU, n, dataInfo, utahMaps] = organize_MEAephysData(datasetpaths,'SU');

%% Params

params.contamthreshold = 10000000000;
params.binsize = 0.025;
params.maxlag = 0.5;
params.fs = 3e4;
params.datasets = 4;

%% Compute the mean firing rates to set a threshold

for i = 1:params.datasets
    
    nUnits = length(SU(i).spikesByTime_condSep);
    
    for units = 1:nUnits
        
        spikenumc = 0;
        
        for cond = 1:8
            
            nTrials = length(SU(i).spikesByTime_condSep{units}.spikesSOAligned{cond});
            
            spikenumtr = 0;
            for tr = 1:nTrials
                
                sn = SU(i).spikesByTime_condSep{units}.spikesSOAligned{cond}{tr}>0;
                spikenumtr = spikenumtr+sum(sn);
            end
            
            spikenumc = spikenumc+(spikenumtr/nTrials);
            
        end
        
        fr(units) = (spikenumc/8)/1;
        
    end
    
    firingRates{i} = fr;
    clear fr
end

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

SUspikesByTimeT = SU(i).spikesByTime_condSep;
SUspikesByTimeIT = SU(i).spikesByTime_interTr_condSep;
cell_pairs{i} = findValPairs3(dist,SU_PWdistances,SUind);
%corrsTrial{i} = sccgparnewcondfr(params,SUspikesByTimeT,sitePairs);
%corrsInterTrial{i} = sccgparnewIT(params,SUspikesByTimeIT,sitePairs);

%corrsTrial{i} = sccgparnewcondfr(params,SUspikesByTimeT,sitePairs);
%corrsInterTrial{i} = sccgparnewcondfrIT(params,SUspikesByTimeIT,sitePairs);

%corrsTrial{i} = vanrossumsync(params,SUspikesByTimeT,sitePairs,'trial');
%corrsInterTrial{i} = vanrossumsync(params,SUspikesByTimeIT,sitePairs,'intertrial');

corrsTrial{i} = vpsync(params,SUspikesByTimeT,sitePairs,'trial');
corrsInterTrial{i} = vpsync(params,SUspikesByTimeIT,sitePairs,'intertrial');


%u = 1:nunits;
%v = (firingRates{i}>=5);
%validFr = u(v);

%corrsTrial{i} = sccgparnewfr(params,SUspikesByTimeT,sitePairs,validFr);
%corrsInterTrial{i} = sccgparnewITfr(params,SUspikesByTimeIT,sitePairs);

%corrsTrial{i} = sccgparnewcondfr(params,SUspikesByTimeT,sitePairs);
end

cd('L:\projects\AbhilashD\Correlations\NewAwakeData')
save('AwakeCorrs_VictorPurpuraDistance.mat','cell_pairs','corrsTrial','corrsInterTrial','datasetpaths','dataInfo','utahMaps');

%% Plot

dBins = 0.5:0.5:4;
lags = linspace(-0.5,0.5,501);

for i = 1:params.datasets
for j = 1:length(cell_pairs{i})
for k = 1:size(cell_pairs{i}{j},1)
pc = corrsInterTrial{i}.meanCCGs{cell_pairs{i}{j}(k,1),cell_pairs{i}{j}(k,2)};
if sum(any(pc==Inf,1))==1 || sum(any(pc==-Inf,1))
    pc = nan(1,501);
end
pairCorrs(k,:) = pc;

end
pairCorrs(pairCorrs==Inf)=0;
pairCorrs(pairCorrs==-Inf)=0;
corrsdBin(j,:) = nanmean(pairCorrs,1);
clear pairCorrs
end
dataSetBinCorrs{i} = corrsdBin;
end

fignames{1} = ['H07_19102016'];
fignames{2} = ['H07_25102016'];
fignames{3} = ['A11_30112016'];
fignames{4} = ['A11_03122016'];

for i = 1:4
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
title([fignames{i} 'Task Synchrony by Distance'])

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
title([fignames{i} 'Task Structure of Synchrony'])

print('-dpng', '-r600', strcat([fignames{i} '5HzCondGM_StructSynchrony_Integral_IT'],'.png'));
end

%% Plot van rossum distances

fignames{1} = ['H07_19102016'];
fignames{2} = ['H07_25102016'];
fignames{3} = ['A11_30112016'];
fignames{4} = ['A11_03122016'];

dBins = 0.5:0.5:4;
lags = linspace(-0.5,0.5,501);

for i = 1:params.datasets
for j = 1:length(cell_pairs{i})
for k = 1:size(cell_pairs{i}{j},1)
pc = corrsTrial{i}.vrsync(cell_pairs{i}{j}(k,1),cell_pairs{i}{j}(k,2));
pairCorrs(k) = pc;

end
corrsdBin(j) = nanmean(pairCorrs);
clear pairCorrs
end
dataSetBinCorrs{i} = corrsdBin;
end

for i = 1:4
figure(i)

plot(dBins,1./dataSetBinCorrs{i})
grid on
axis tight
box off
xlabel('distance in mm')
ylabel('1/VPD - Synchrony')
title([fignames{i} 'Task Structure of Synchrony'])

print('-dpng', '-r600', strcat([fignames{i} 'inv_VPD_StructSynchrony'],'.png'));
end