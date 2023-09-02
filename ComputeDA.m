clear all
clc
%parpool(8)
matlabpool open 8

%% Load datasets

load('spikeData_WithPWstuff_20150708')

%% Set parameters

params.fulltime = 5000; % in samples for on and off
params.maxlag = 5000; % for Sccg
params.numtrials=200;
params.intwin = 500; % for computing r_ccg
tag = [1 0]; % Use tag(ind) in the function argument to compute the CCGS for either ON or OFF phases
params.contamthreshold = 0.1; % increase or decrease this to change the number of WI units
params.datsets = 4; % number of datasets
params.firingthreshold = 1; % Mathematically, one needs 1 spike in a timeseries to compute a crosscorrelogram
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
    
    spikeTime = SU(1,i).spikesBySampleIndex;
    
    spikeTimes = cell(length(spikeTime),params.numtrials);
    
    for kk = 1:length(spikeTime)
        for j = 1:params.numtrials
            spikeTimes{kk,j}=spikeTime{kk}{j}; % spike times
        end
    end
    
    SUcontam = utahMaps(1,i).SU.integratedTabel(:,5);
    
    SU_angles = SU(1,i).PWanglesFromXaxis;
    
    SU_PWdistances = SU(1,i).PWdistances;
    
    indmat = ones(size(SU_angles));
    indmat = tril(indmat);
    SU_PWdistances = SU_PWdistances(indmat==1);
    SU_PWdistances(indmat==0)=nan;% Remove doubles
    
    SU_angles = SU_angles(indmat==1);
    SU_angles(indmat==0)=nan;% Remove doubles
    
    distangs = identifypairsangs(SU_angles,SU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    SUind = find(SUcontam<=params.contamthreshold);
    
    % Phase? Stim on or stim off_
    
    phase = tag(2); % On
    
    % Eliminate low spiking cells
    
    cell_pairs_angs = findValPairs(params,distangs,spikeTimes,SUind,phase);
    
    % Start computations
    
    %Compute shuffle-corrected cross-correlograms
    
    disp('Calculation of SC starts')
    
    for u = 1:size(distangs,1)
        parfor v = 1:size(distangs,2)
            if isempty(cell_pairs_angs{u,v})==1
                rccgangs{u,v}=0;
            else
                [rccgangs{u,v}] = sccgpar_distang(params,spikeTimes,cell_pairs_angs{u,v},phase);
            end
        end
    end
    

    resultsSCangs{i}.Rccg = rccgangs;
    disp('calculation of SC is done')
    
     disp('Calculation of SSC starts')
    for u = 1:size(distangs,1)
        parfor v = 1:size(distangs,2)
            if isempty(cell_pairs_angs{u,v})==1
                sscangs{u,v}=0;
            else
                [sscangs{u,v}, freqs] = spikespikecoherence_distang(params,cell_pairs_angs{u,v},spikeTimes,phase);
            end
        end
    end
    
    resultsSSCangs{i}=sscangs;
    
    disp('Calculation of SSC is done')
    clc
    msg=sprintf('Data set %d of %d is done...\n',i,nDataSet);
    disp(msg);
    
    endTime=toc;
    msg=sprintf('Total elapsed time for %dth dataset is %d seconds',i,endTime);
    disp(msg)
    
end

save('angsSUSCallOff3GMp01.mat','resultsSCangs', '-v7.3')
save('angsSUallOffSSC3p0110150.mat','resultsSSCangs', '-v7.3')

%delete(gcp)

matlabpool close
