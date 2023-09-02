contamind = [];

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
    
    dist{i} = identifypairs(SU_PWdistances);
    
    % Identify well-isolated units by the indices of contamination
    
    SUind = find(SUcontam<=params.contamthreshold);
    contamind = [contamind; SUcontam(SUind)];
    
    %cell_pairs{i} = findValPairs3(dist{i},SU_PWdistances,SUind);
    
end