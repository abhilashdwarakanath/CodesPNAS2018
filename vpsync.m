function corrs = vpsync(params,SUspikesByTime,pairs,tag)

tf = 500;
decimFac = tf/params.fs;
T = 1/tf:1/tf:1.3;

% Do first for 90 to 270

for i = 1:size(pairs,1)
    cell_pair = pairs(i,:);
    fprintf('Starting a cell pair...\n');
    
    switch tag
        
        case 'trial'
            
            cell1 = SUspikesByTime{cell_pair(1)}.spikesSOAligned;
            cell2 = SUspikesByTime{cell_pair(2)}.spikesSOAligned;
            
        case 'intertrial'
            cell1 = SUspikesByTime{cell_pair(1)}.spikesEOAligned;
            cell2 = SUspikesByTime{cell_pair(2)}.spikesEOAligned;
    end
    tic;
    for k = 1:length(cell1)
        fprintf('Computing correlations for condition %d\n',k);
        
        for l=1:length(cell1{k})
            fprintf('Trial %d of %d\n',l,length(cell1{k}));
            spktimes1 = ceil((cell1{k}{l}.*30)*decimFac);
            spktimes1 = spktimes1(spktimes1>0);
            spktimes2 = ceil((cell2{k}{l}.*30)*decimFac);
            spktimes2 = spktimes2(spktimes2>0);
            
            if ~isempty(spktimes1) && ~isempty(spktimes2)
            
                cc1(l) = compute_normalized_dist(spktimes1,spktimes2,2);
            
            else
                cc1(l) = NaN;
            end
            
        end
        toc;
        
        synchrony{cell_pair(1),cell_pair(2)}{k} = nanmean(cc1);
        
    end
    
    for k = 1:length(cell1)
        tcc(k) = synchrony{cell_pair(1),cell_pair(2)}{k};
        
    end
    
    mvrsync(cell_pair(1),cell_pair(2)) = nanmean(tcc);
    
    
end

corrs.vrsync = mvrsync;