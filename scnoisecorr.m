%clear all
close all
clc

%% Load data

load('spikedata_151118.mat')

%% Set parameters

params.fulltime = 5000; % in samples for on and off
params.maxlag = 2000; % for Sccg
params.numtrials=200;
params.jmaxlag = 2000; % for Jccg
params.intwin = [5 25 50 250 500 1000 2000]; % for computing r_ccg
tag = [1 0]; % Use tag(ind) in the function argument to compute the CCGS for either ON or OFF phases
params.contamthreshold = 0.05; % increase or decrease this to change the number of WI units
params.jitwin = [5 10 25 50];% in samples. jitwin*2 = ms
params.datsets = 4; % number of datasets
params.firingthreshold = 0; % Mathematically, one needs 1 spike in a timeseries to compute a crosscorrelogram
params.fpass = [1 100];
params.trialave = 1;
params.Fs = 500;
params.tapers = [];
params.err = 1;
movingwin = [200e-3 20e-3];


% This program is a rough code to compute shuffle-corrected noise
% correlations (which doesn't really make sense because one always removes
% the mean from the trace and then computes noise corrs)

%% Do Single Units
% for z = 1:4
%     [SUncallon{z}, SUncalloff{z}] = shufflednoisecorr(SU,z);
% end

for z = 1:4
    nDataSet = z;
    msg=sprintf('Data set %d of %d starts...\n',z,nDataSet);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances and contamination values
    
    spikeTime = SU(1,z).spikesBySampleIndex;
    
    spikeTimes = cell(length(spikeTime),params.numtrials);
    
    for kk = 1:length(spikeTime)
        for j = 1:params.numtrials
            spikeTimes{kk,j}=spikeTime{kk}{j}; % spike times
        end
    end
    
    SU_PWdistances = SU(1,z).PWdistances; % distances
    
    SUcontam = utahMaps(1,z).SU.integratedTabel(:,5); % contamination indices
    
    SU_PWdistances = tril(SU_PWdistances); % Remove doubles
    
    % Identify pairs in distance bins
    
    dist = identifypairs(SU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    SUind = find(SUcontam<=params.contamthreshold);
    %SUind = find(SUcontam<=100000000);
    
    % Eliminate low spiking cells
    
    [SUrawncallon{z}, SUrawncalloff{z}] = noisecorr(SU,z);
    
    SUcell_pairs{z} = findValPairs3(dist,SUrawncallon{z},SUind);
end
% 
% for z = 1:4
%     [frON{z},frOFF{z}] = psth(SU,z);
% end

%% Do MUA
% for z = 1:4
%     [MUncallon{z}, MUncalloff{z}] = shufflednoisecorr(SU,z);
% end

for z = 1:4
     nDataSet = z;
    msg=sprintf('Data set %d of %d starts...\n',z,nDataSet);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances and contamination values
    
    spikeTime = MU(1,z).spikesBySampleIndex;
    
    spikeTimes = cell(length(spikeTime),params.numtrials);
    
    for kk = 1:length(spikeTime)
        for j = 1:params.numtrials
            spikeTimes{kk,j}=spikeTime{kk}{j}; % spike times
        end
    end
    
    MU_PWdistances = MU(1,z).PWdistances; % distances
    
    MUcontam = utahMaps(1,z).MU.integratedTabel(:,5); % contamination indices
    
    MU_PWdistances = tril(MU_PWdistances); % Remove doubles
    
    % Identify pairs in distance bins
    
    dist = identifypairs(MU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    MUind = find(MUcontam<=params.contamthreshold);
    
    % Eliminate low spiking cells
    
    [MUrawncallon{z}, MUrawncalloff{z}] = noisecorr(MU,z);
    
    MUcell_pairs{z} = findValPairs3(dist,MU_PWdistances,MUind);
end

% for z = 1:4
%     [frON{z},frOFF{z}] = psth(SU,z);
% end
%% Do jMU

% for z = 1:4
%     [ncallon{z}, ncalloff{z}] = shufflednoisecorr(SU,z);
% end

for z = 1:4
     nDataSet = z;
    msg=sprintf('Data set %d of %d starts...\n',z,nDataSet);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances and contamination values
    
    spikeTime = jMU(1,z).spikesBySampleIndex;
    
    spikeTimes = cell(length(spikeTime),params.numtrials);
    
    for kk = 1:length(spikeTime)
        for j = 1:params.numtrials
            spikeTimes{kk,j}=spikeTime{kk}{j}; % spike times
        end
    end
    
    jMU_PWdistances = MU(1,z).PWdistances; % distances
    
    jMUcontam = utahMaps(1,z).MU.integratedTabel(:,5); % contamination indices
    
    jMU_PWdistances = tril(jMU_PWdistances); % Remove doubles
    
    % Identify pairs in distance bins
    
    dist = identifypairs(jMU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    jMUind = find(jMUcontam<=params.contamthreshold);
    
    % Eliminate low spiking cells
    
    [jMUrawncallon{z}, jMUrawncalloff{z}] = noisecorr(jMU,z);
    
    jMUcell_pairs{z} = findValPairs3(dist,jMU_PWdistances,jMUind);
end

%% Compute correlation structure
dBins = 0.5:0.5:4;
% SU

[SUnc] = corrstruct3(params,SUcell_pairs,dBins,SUrawncallon,SUrawncallon);

% MU
[MUnc] = corrstruct3(params,MUcell_pairs,dBins,MUrawncallon,MUrawncallon);

% jMU
[jMUnc] = corrstruct3(params,jMUcell_pairs,dBins,jMUrawncallon,jMUrawncallon);
