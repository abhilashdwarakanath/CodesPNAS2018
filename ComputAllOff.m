clear all
clc
matlabpool open 8
%% Load datasets
load('spikeData20150618_UATE2010.mat')

%% Generate shuffled correlation
resultsSC = cell(1,4);
resultsJC = cell(1,4);
x = cell(1,4);

for y = 1:4
    resultsJC{y} = x;
end

%jitwin=[10 20 50 100];
jitwin=[5 10 25 50];% in ms. jitwin/2 = ms

for i = 1:4
    
    nDataSet = i;
    msg=sprintf('Data set %d of %d starts...\n',i,nDataSet);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances and contamination values
    
    spikeTime = SU(1,i).spikesBySampleIndex;
    
    for kk = 1:length(spikeTime)
        for j = 1:200
            spikeTimes{kk,j}=spikeTime{kk}{j};
        end
    end
    
    SU_PWdistances = SU(1,i).PWdistances;
    SUcontam = utahMaps(1,i).SU.integratedTabel(:,5);
    
    SU_PWdistances = tril(SU_PWdistances);
    
    % Identify pairs in distance bins
    
    dist = identifypairs(SU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    SUind = find(SUcontam<=0.05);
    
    % Eliminate low spiking cells
    
    cell_pairs = findValPairs(dist,spikeTimes,SUind);
    
    % Start computations
   
    disp('Calculation of SC starts')
    parfor z = 1:length(dist)
        [ccSC{z}, rccg{z} A1{z} A2{z}] = sccgparoff(spikeTimes,cell_pairs,z);
    end
    
    resultsSC{i}.SC = ccSC;
    resultsSC{i}.Rccg = rccg;
    resultsSC{i}.A1 = A1;
    resultsSC{i}.A2 = A2;
    resultsSC{i}.Rccg = rccg;
    disp('calculation of SC is done')
    
    for tt = 1:length(jitwin)
        disp('Calculation of JC starts...')
        parfor z = 1:length(dist)
            [ccJC{z}] = jccgparoff(spikeTimes,cell_pairs,jitwin,z)
        end
        resultsJC{i}{tt}.JC = ccJC;
        disp('Calculation of JC is done...')
    end
    
    clc
    msg=sprintf('Data set %d of %d is done...\n',i,nDataSet);
    disp(msg);
    
    endTime=toc;
    msg=sprintf('Total elapsed time for %dth dataset is %d seconds',i,endTime);
    disp(msg)
    
end

matlabpool close
