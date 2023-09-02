function [dataset,rscON,rscOFF] = corrstruct2(dist,params,dBins,aBins,Qon,Qoff,dataset)

rscON = zeros(length(dBins),length(aBins));
rscOFF = zeros(length(dBins),length(aBins));

for currset = 1:params.ndatasets
    
    % Compute correlation coefficients along distance
    
    corr_dist_on = zeros(size(dBins,2),size(aBins,2));
    corr_dist_off = zeros(size(dBins,2),size(aBins,2));
    
    %varon = diag(sqrt(diag((Qon{currset}))));
    %varoff = diag(sqrt(diag((Qoff{currset}))));
    
    corrmon = (Qon{currset});%inv(varon)*Qon{currset}*inv(varon);
    corrmoff = (Qoff{currset});%inv(varoff)*Qoff{currset}*inv(varoff);
    
    distsemon = zeros(length(dBins),length(aBins));
    distsemoff = zeros(length(dBins),length(aBins));
    
    for nBins = 1:size(dBins,2)
        for nAngs = 1:size(aBins,2)
            erron = [];
            erroff = [];
            for nPairs = 1:length(dist{currset}{nBins,nAngs})
                
                if ~isempty(dist{currset}{nBins,nAngs})
                    
                    [x,y] = ind2sub(size(Qon{currset}),dist{currset}{nBins,nAngs}(nPairs));
                    erron=[erron corrmon(x,y)];
                    erroff=[erroff corrmoff(x,y)];
                    corr_dist_on(nBins,nAngs) = corr_dist_on(nBins,nAngs)+corrmon(x,y);
                    corr_dist_off(nBins,nAngs) = corr_dist_off(nBins,nAngs)+corrmoff(x,y);
                    
                else
                    erron=[erron nan];
                    erroff=[erroff nan];
                    corr_dist_on(nBins,nAngs) = nan;
                    corr_dist_off(nBins,nAngs) = nan;
                    
                end
                
            end
            
            distsemon(nBins,nAngs) = nanstd(erron)/sqrt(length(erron));
            distsemoff(nBins,nAngs) = nanstd(erroff)/sqrt(length(erroff));
            
            if ~isempty(dist{currset}{nBins,nAngs})
                corr_dist_on(nBins,nAngs) = corr_dist_on(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
                corr_dist_off(nBins,nAngs) = corr_dist_off(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
                
                distsemon(nBins,nAngs) = distsemon(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
                distsemoff(nBins,nAngs) = distsemoff(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
            else
                corr_dist_on(nBins,nAngs) = nan;
                corr_dist_off(nBins,nAngs) = nan;
                
                distsemon(nBins,nAngs) = nan;
                distsemoff(nBins,nAngs) = nan;
            end
            
        end
        
    end
    
    if length(dBins) > 1 && length(aBins) > 1
        dataset(currset).modelON.nc = corr_dist_on;
        dataset(currset).modelOFF.nc = corr_dist_off;
        dataset(currset).modelON.ncerr = distsemon;
        dataset(currset).modelOFF.ncerr = distsemoff;
    elseif length(dBins) > 1 && length(aBins) == 1
        dataset(currset).modelON.ncdist = corr_dist_on;
        dataset(currset).modelOFF.ncdist = corr_dist_off;
        dataset(currset).modelON.ncdisterr = distsemon;
        dataset(currset).modelOFF.ncdisterr = distsemoff;
    elseif length(dBins) == 1 && length(aBins) > 1
        dataset(currset).modelON.ncangs = corr_dist_on;
        dataset(currset).modelOFF.ncangs = corr_dist_off;
        dataset(currset).modelON.ncangserr = distsemon;
        dataset(currset).modelOFF.ncangserr = distsemoff;
    else
        dataset(currset) = dataset(currset);
    end
    
    rscON(:,:,currset) = corr_dist_on;
    rscOFF(:,:,currset) = corr_dist_off;
    
end

rscON = nanmean(rscON,3);
rscOFF = nanmean(rscOFF,3);
