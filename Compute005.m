clear all
clc

%% Load datasets

load('spikedata_151118.mat')

%% Set parameters

params.fulltime = 5000; % in samples for on and off
params.maxlag = 5000; % for Sccg
params.numtrials=200;
params.jmaxlag = 1000; % for Jccg
params.intwin = [5 50 250 500 2500 5000]; % for computing r_ccg
tag = [1 0]; % Use tag(ind) in the function argument to compute the CCGS for either ON or OFF phases
params.contamthreshold = 100000; % increase or decrease this to change the number of WI units
params.jitwin = [5 10 25 50];% in samples. jitwin*2 = ms
params.datsets = 4; % number of datasets
params.firingthreshold = 1; % Mathematically, one needs 1 spike in a timeseries to compute a crosscorrelogram
params.fpass = [0.1 250];
params.trialave = 1;
params.Fs = 500;
params.tapers = [5 2];
params.err = 1;
movingwin = [200e-3 20e-3];
%% Do computations and collect the results

ONresultsSC = cell(1,params.datsets);
%ONresultsSJC = cell(1,params.datsets);
%ONresultsSSC= cell(1,params.datsets);
%ONresultsJC = cell(1,params.datsets);

OFFresultsSC = cell(1,params.datsets);
%OFFresultsSJC = cell(1,params.datsets);
%OFFresultsSSC= cell(1,params.datsets);
%OFFresultsJC = cell(1,params.datsets);

x = cell(1,params.datsets);

% for y = 1:params.datsets
%     ONresultsJC{y} = x;
%     OFFresultsJC{y} = x;
% end

for i = 1:params.datsets
    
    nDataSet = i;
    msg=sprintf('Data set %d of %d starts...\n',i,nDataSet);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances and contamination values
    
    spikeTime = SU(1,i).spikesBySampleIndex;
    
    spikeTimes = cell(length(spikeTime),params.numtrials);
    
    for kk = 1:length(spikeTime)
        for j = 1:params.numtrials
            spikeTimes{kk,j}=spikeTime{kk}{j}; % spike times
        end
    end
    
    SU_PWdistances = SU(1,i).PWdistances; % distances
    
    SUcontam = utahMaps(1,i).SU.integratedTabel(:,5); % contamination indices
    
    SU_PWdistances = tril(SU_PWdistances); % Remove doubles
    
    % Identify pairs in distance bins
    
    dist = identifypairs(SU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    SUind = find(SUcontam<=params.contamthreshold);
    
    % ON PHASE
    
    phase = tag(1); % on
    
    % Eliminate low spiking cells
    
    cell_pairs = findValPairs(params,dist,spikeTimes,SUind,phase);
    
    % Start computations
    
    %Compute shuffle-corrected cross-correlograms
    disp('Calculation of SC starts')
    parfor z = 1:length(dist)
        [ccSC(:,z), rccg(:,z),rccg_tb(:,z)] = sccgpar(params,spikeTimes,cell_pairs,z,phase);
    end
    
    ONresultsSC{i}.SC = ccSC;
    ONresultsSC{i}.Rccg = rccg;
    ONresultsSC{i}.Rccg_tb = rccg_tb;
    disp('calculation of SC is done')
    
    % Compute Shuffle and Jitter corrected CCGs
    
%     parfor z = 1:length(dist)
%         [ccSJC(:,z)] = sjccgpar(params,spikeTimes,cell_pairs,params.jitwin(3),z,phase);
%     end
%     
%     ONresultsSJC{i}.SJC = ccSJC;
    
    % compute jitter-corrected cross-correlograms
    
%     for tt = 1:length(params.jitwin)
%         disp('Calculation of JC starts...')
%         jw = params.jitwin(tt);
%         parfor z = 1:length(dist)
%             [ccJC(:,z)] = jccgpar(params,spikeTimes,cell_pairs,jw,z,phase)
%         end
%         ONresultsJC{i}{tt}.JC = ccJC;
%         disp('Calculation of JC is done...')
%     end
%     
    parfor z = 1:length(dist)
        [ssc(:,z), freqs] = spikespikecoherence(params,cell_pairs,spikeTimes,z,phase);
    end
    
    ONresultsSSC{i}=ssc;
    
    clc
    msg=sprintf('Data set %d of %d is done...\n',i,nDataSet);
    disp(msg);
    
    endTime=toc;
    msg=sprintf('Total elapsed time for %dth dataset is %d seconds',i,endTime);
    disp(msg)
    
    % OFF PHASE
    
    phase = tag(2); % off
    
    % Eliminate low spiking cells
    
    cell_pairs = findValPairs(params,dist,spikeTimes,SUind,phase);
    
    % Start computations
    
    %Compute shuffle-corrected cross-correlograms
    disp('Calculation of SC starts')
    parfor z = 1:length(dist)
        [ccSC(:,z), rccg(:,z),rccg_tb(:,z)] = sccgpar(params,spikeTimes,cell_pairs,z,phase);
    end
    
    OFFresultsSC{i}.SC = ccSC;
    OFFresultsSC{i}.Rccg = rccg;
    OFFresultsSC{i}.Rccg_tb = rccg_tb;
    disp('calculation of SC is done')
    
    % Compute Shuffle and Jitter corrected CCGs
    
%     parfor z = 1:length(dist)
%         [ccSJC(:,z)] = sjccgpar(params,spikeTimes,cell_pairs,params.jitwin(3),z,phase);
%     end
%     
%     OFFresultsSJC{i}.SJC = ccSJC;
%     
%     % compute jitter-corrected cross-correlograms
%     
%     for tt = 1:length(params.jitwin)
%         disp('Calculation of JC starts...')
%         jw = params.jitwin(tt);
%         parfor z = 1:length(dist)
%             [ccJC(:,z)] = jccgpar(params,spikeTimes,cell_pairs,jw,z,phase)
%         end
%         OFFresultsJC{i}{tt}.JC = ccJC;
%         disp('Calculation of JC is done...')
%     end
%     
    parfor z = 1:length(dist)
        [ssc(:,z), freqs] = spikespikecoherence(params,cell_pairs,spikeTimes,z,phase);
    end
    
    OFFresultsSSC{i}=ssc;
    
    clc
    msg=sprintf('Data set %d of %d is done...\n',i,nDataSet);
    disp(msg);
    
    endTime=toc;
    msg=sprintf('Total elapsed time for %dth dataset is %d seconds',i,endTime);
    disp(msg)
    
end

resultsON.shufflecorr = ONresultsSC;
%resultsON.sjcorr = ONresultsSJC;
%resultsON.jittercorr = ONresultsJC;
%resultsON.spikecoher = ONresultsSSC;

resultsOFF.shufflecorr = OFFresultsSC;
%resultsOFF.sjcorr = OFFresultsSJC;
%resultsOFF.jittercorr = OFFresultsJC;
%resultsOFF.spikecoher = OFFresultsSSC;

delete(gcp)
