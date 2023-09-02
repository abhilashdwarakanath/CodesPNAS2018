function [nc] = corrstruct3(params,dist,dBins,Qon,Qoff)

for currset = 1:params.ndatasets
    
    % Compute correlation coefficients along distance
    
    corr_dist_on = zeros(1,size(dBins,2));
    corr_dist_off = zeros(1,size(dBins,2));
    
    %corrmon = corrcov((Qon{currset}));
    %corrmoff = corrcov((Qoff{currset}));
    
    corrmon = Qon{currset};
    corrmoff = Qoff{currset};
    
    distsemon = zeros(1,length(dBins));
    distsemoff = zeros(1,length(dBins));
    
    
    for nBins = 1:size(dBins,2)
        erron = [];
        erroff = [];
        for nPairs = 1:length(dist{currset}{nBins})
            
            
            [inds(1),inds(2)] = ind2sub(size(Qon{currset}),dist{currset}{nBins}(nPairs));
            %inds=dist{currset}{nBins}(nPairs,:);
            erron=[erron corrmon(inds(1),inds(2))];
            erroff=[erroff corrmoff(inds(1),inds(2))];
            
        end
        
        distsemon(nBins) = nanstd(erron)/sqrt(length(erron));
        distsemoff(nBins) = nanstd(erroff)/sqrt(length(erroff));
        
        corr_dist_on(nBins) = nanmean(erron);
        corr_dist_off(nBins) = nanmean(erroff);
        
    end

       nc.ON(:,currset) = corr_dist_on;
       nc.OFF(:,currset) = corr_dist_off;
       nc.errON(:,currset) = distsemon;
       nc.errOFF(:,currset) = distsemoff;
    
end

