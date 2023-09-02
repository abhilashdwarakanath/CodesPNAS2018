clear all
clc
parpool(8)

%% Load datasets

load('spikeData_WithPWstuff_20150708')

%% Set parameters

params.fulltime = 5000; % in samples for on and off
params.maxlag = 5000; % for Sccg
params.numtrials=200;
params.intwin = [1 5 10 20 50 100 125 250 500]; % for computing r_ccg
tag = [1 0]; % Use tag(ind) in the function argument to compute the CCGS for either ON or OFF phases
params.datsets = 4; % number of datasets
params.firingthreshold = 3; % Mathematically, one needs 1 spike in a timeseries to compute a crosscorrelogram
params.fpass = [10 150];
params.trialave = 1;
params.Fs = 500;
params.tapers = [];
params.err = 0;
movingwin = [200e-3 20e-3];
%% Do computations and collect the results

resultsSCangs = cell(1,params.datsets);
resultsSSCangs = cell(1,params.datsets);

for i = 1:params.datsets
    
    nDataSet = i;
    msg=sprintf('Data set %d of %d starts...\n',i,nDataSet);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances and contamination values
    
    spikeTime = MU(1,i).spikesBySampleIndex;
    
    spikeTimes = cell(length(spikeTime),params.numtrials);
    
    for kk = 1:length(spikeTime)
        for j = 1:params.numtrials
            spikeTimes{kk,j}=spikeTime{kk}{j}; % spike times
        end
    end
    
    MU_angles = MU(1,i).PWanglesFromXaxis;
    
    MU_angles = tril(MU_angles);
    
    % Identify pairs in distance bins

    cell_pairs_angs = identifypairsangs(MU_angles);
    
    % Phase? Stim on or stim off_
    
    phase = tag(2); % On
    
    % Start computations
    
    %Compute shuffle-corrected cross-correlograms
    
        disp('Calculation of SC starts')

        parfor z = 1:length(cell_pairs_angs)
            [ccSCangs{z}, rccgangs{z}] = sccgpar(params,spikeTimes,cell_pairs_angs,z,phase);
        end
    
        resultsSCangs{i}.SC = ccSCangs;
        resultsSCangs{i}.Rccg = rccgangs;
        disp('calculation of SC is done')
    
    parfor z = 1:length(cell_pairs_angs)
        [sscangs{z}, freqs] = spikespikecoherence(params,cell_pairs_angs,spikeTimes,z,phase);
    end
    
    resultsSSCangs{i}=sscangs;
    
    clc
    msg=sprintf('Data set %d of %d is done...\n',i,nDataSet);
    disp(msg);
    
    endTime=toc;
    msg=sprintf('Total elapsed time for %dth dataset is %d seconds',i,endTime);
    disp(msg)
    
end

save('angsMUSCallOff3GMp01.mat','resultsSCangs', '-v7.3')
save('angsMUallOffSSC3p0110150.mat','resultsSSCangs', '-v7.3')

delete(gcp)
