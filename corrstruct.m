function [dataset,corrstrucdistON,corrstrucdistOFF,rem_corrstrucdistON,rem_corrstrucdistOFF] = corrstruct(dist,params,dBins,aBins,dataset)

for currset = 1:params.ndatasets
    
    residuals_on = (dataset(currset).modelON.R);
    residuals_off = (dataset(currset).modelOFF.R);
    rem_on = corrcov(dataset(currset).modelON.rem);
    rem_off = corrcov(dataset(currset).modelOFF.rem);
    
    C_on =  (dataset(currset).modelON.C);
    C_on =C_on*C_on';
    C_on = (C_on+residuals_on);
    
    C_off =  (dataset(currset).modelOFF.C);
    C_off =C_off*C_off';
    C_off = (C_off+residuals_off);
    %varon = diag(sqrt(diag((C_on))));
    %varoff = diag(sqrt(diag((C_off))));
    
    %Corr_on=inv(varon)*C_on*inv(varon);
    %Corr_off=inv(varoff)*C_off*inv(varoff);
    Corr_on = corrcov(C_on);%inv(varon)*C_on*inv(varon);
    Corr_off = corrcov(C_off); %inv(varoff)*C_off*inv(varoff);
    
    distax = -1:0.05:1;
    lda = length(distax);
    
    if length(dBins)>2 && length(aBins)>2
        
        figure(100+currset)

        subplot(2,2,3)
        histogram(rem_on(:),distax)
        title('Stim ON Residual')
        
        subplot(2,2,4)
        histogram(rem_off(:),distax)
        title('Stim Off Residual')
        
        figure(200+currset)
        subplot(2,2,1)
        imagescnan(Corr_on)
        title('Stim ON LV')
        
        subplot(2,2,2)
        imagescnan(Corr_off)
        title('Stim Off LV')
        
        subplot(2,2,3)
        imagescnan(rem_on)
        title('Stim ON Residual')
        
        subplot(2,2,4)
        imagescnan(rem_off)
        title('Stim Off Residual')
        
    end
    
    
    % Compute correlation coefficients along distance
    
    distsemon = zeros(length(dBins),length(aBins));
    distsemoff = zeros(length(dBins),length(aBins));
    
    corr_dist_on = zeros(size(dBins,2),size(aBins,2));
    corr_dist_off = zeros(size(dBins,2),size(aBins,2));
    rem_corr_dist_on = zeros(size(dBins,2),size(aBins,2));
    rem_corr_dist_off = zeros(size(dBins,2),size(aBins,2));
    
    for nBins = 1:size(dBins,2)
        for nAngs = 1:size(aBins,2)
            erron = [];
            erroff = [];
            for nPairs = 1:length(dist{currset}{nBins,nAngs})
                
                if ~isempty(dist{currset}{nBins,nAngs})
                    
                    [x,y] = ind2sub(size(C_on),dist{currset}{nBins,nAngs}(nPairs));
                    erron=[erron Corr_on(x,y)];
                    erroff=[erroff Corr_off(x,y)];
                    corr_dist_on(nBins,nAngs) = corr_dist_on(nBins,nAngs)+Corr_on(x,y);
                    corr_dist_off(nBins,nAngs) = corr_dist_off(nBins,nAngs)+Corr_off(x,y);
                    rem_corr_dist_on(nBins,nAngs) = rem_corr_dist_on(nBins,nAngs)+rem_on(x,y);
                    rem_corr_dist_off(nBins,nAngs) = rem_corr_dist_off(nBins,nAngs)+rem_off(x,y);
                else
                    erron=[erron nan];
                    erroff=[erroff nan];
                    corr_dist_on(nBins,nAngs) = nan;
                    corr_dist_off(nBins,nAngs) = nan;
                    rem_corr_dist_on(nBins,nAngs) = nan;
                    rem_corr_dist_off(nBins,nAngs) = nan;
                    
                end
                
            end
            
            distsemon(nBins,nAngs) = nanstd(erron)/sqrt(length(erron));
            distsemoff(nBins,nAngs) = nanstd(erroff)/sqrt(length(erroff));
            
            if ~isempty(dist{currset}{nBins,nAngs})
                rem_corr_dist_on(nBins,nAngs) = rem_corr_dist_on(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
                rem_corr_dist_off(nBins,nAngs) = rem_corr_dist_off(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
                corr_dist_on(nBins,nAngs) = corr_dist_on(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
                corr_dist_off(nBins,nAngs) = corr_dist_off(nBins,nAngs)/length(dist{currset}{nBins,nAngs});
            else
                rem_corr_dist_on(nBins,nAngs) = nan;
                rem_corr_dist_off(nBins,nAngs) = nan;
                corr_dist_on(nBins,nAngs) = nan;
                corr_dist_off(nBins,nAngs) = nan;
                
            end
            
        end
        
    end
    
    if length(dBins) > 2 && length(aBins) > 2
        dataset(currset).modelON.corrstruc = corr_dist_on;
        dataset(currset).modelOFF.corrstruc = corr_dist_off;
        dataset(currset).modelON.rem_corrstruc = rem_corr_dist_on;
        dataset(currset).modelOFF.rem_corrstruc = rem_corr_dist_off;
        dataset(currset).modelON.sem = distsemon;
        dataset(currset).modelOFF.sem = distsemoff;
        
    else
        dataset(currset) = dataset(currset);
    end
    
    corrstrucdistON(:,:,currset) = corr_dist_on;
    corrstrucdistOFF(:,:,currset) = corr_dist_off;
    
    semON(:,:,currset) = distsemon;
    semOFF(:,:,currset) = distsemoff;
    
    rem_corrstrucdistON(:,:,currset) = rem_corr_dist_on;
    rem_corrstrucdistOFF(:,:,currset) = rem_corr_dist_off;
    
    
end

corrstrucdistON = nanmean(corrstrucdistON,3);
corrstrucdistOFF = nanmean(corrstrucdistOFF,3);
rem_corrstrucdistON = nanmean(rem_corrstrucdistON,3);
rem_corrstrucdistOFF = nanmean(rem_corrstrucdistOFF,3);