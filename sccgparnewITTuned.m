function corrs = sccgparnewITTuned(params,SUspikesByTime,pairs,pvals)

tf = 500;
decimFac = tf/params.fs;
T = 1/tf:1/tf:1.3;

% Do first for 90 to 270

for i = 1:size(pairs,1)
    cell_pair = pairs(i,:);
    fprintf('Starting a cell pair...\n');
    cell1 = SUspikesByTime{cell_pair(1)}.spikesEOAligned;
    cell2 = SUspikesByTime{cell_pair(2)}.spikesEOAligned;
    
    p1 = pvals(cell_pair(1));
    p2 = pvals(cell_pair(2));
    
    if p1<=0.05 && p2<=0.05
        
        
    for k = 1:length(cell1)
        fprintf('Computing correlations for condition %d\n',k);
        N = length(cell1{k});
        original_indices = 1:N;
        done = 0;
        while ~done
            shuffling = randperm(N); % this will hold the result when loop is exited
            done = all(shuffling~=original_indices);
        end
        clear corr_cc1; clear peakShiftstt;
        tic;
        for l=1:length(cell1{k})
                        fprintf('Trial %d of %d\n',l,length(cell1{k}));
            spktimes1 = ceil((cell1{k}{l}.*30)*decimFac);
            spktimes1(spktimes1==0) = 1;
            spktimes2 = ceil((cell2{k}{shuffling(l)}.*30)*decimFac);
            spktimes2(spktimes2==0) = 1;
            spktimes3 = ceil((cell2{k}{l}.*30)*decimFac);
            spktimes3(spktimes3==0) = 1;
            
            if numel(spktimes1)>=1 && numel(spktimes2)>=1 && numel(spktimes3)>=1
            
            st1 = zeros(1,length(T));
            st2 = zeros(1,length(T));
            st3 = zeros(1,length(T));
            st1(spktimes1) = 1;
            st2(spktimes2) = 1;
            st3(spktimes3) = 1;
            
            fr1 = sum(st1)/(length(st1)/tf);
            fr2 = sum(st3)/(length(st3)/tf);
            fr3 = sum(st2)/(length(st2)/tf);
            
            [cc1,lags] = cross_corr(st1,st2,params.maxlag,(1/tf)); % Shuffled CCG
            theta = length(lags)./(length(lags)-abs(lags/tf));
            cc1 = cc1./geomean([fr1 fr2]);
            cc1 = cc1./theta;
            [cc2,lags] = cross_corr(st1,st3,params.maxlag,(1/tf)); % Expected CCG
            cc2 = cc2./geomean([fr1 fr3]);
            cc2 = cc2./theta;
            differ = cc2-cc1; %#ok<*PFOUS>
            corr_cc1(l,:) = differ;
            gm(l) =  geomean([fr1 fr3]);
            [~,idx] = max(abs(corr_cc1(l,:)));
            if idx==1                 timeDiff = [0];             else             timeDiff = lags(idx);             end
            peakShiftstt(l) = timeDiff;
            
            else
                
                corr_cc1(l,:) = nan(1,501);
                gm(l) = nan;
                peakShiftstt(l) = nan;
            end
            
        end
        toc;
        
        zeroRows = all(corr_cc1==0,2);
        corr_cc1(zeroRows,:) = [];
        peakShiftstt(zeroRows) = [];
        gm(zeroRows) = [];
        geoMean{cell_pair(1),cell_pair(2)}{k} = nanmean(gm);
        crosscorrs{cell_pair(1),cell_pair(2)}{k} = nanmean(corr_cc1,1);
        peakShifts{cell_pair(1),cell_pair(2)}{k} = nanmean(peakShiftstt);
        
    end
    
    for k = 1:length(cell1)
        tcc(k,:) = crosscorrs{cell_pair(1),cell_pair(2)}{k};
        tgm(k) = geoMean{cell_pair(1),cell_pair(2)}{k};
        tps(k) = peakShifts{cell_pair(1),cell_pair(2)}{k};
    end
    
    mcrosscorrs{cell_pair(1),cell_pair(2)} = nanmean(tcc,1);
    mgm{cell_pair(1),cell_pair(2)} = nanmean(tgm);
    mps{cell_pair(1),cell_pair(2)} = nanmean(tps,1);
    
     else
        mcrosscorrs{cell_pair(1),cell_pair(2)} = NaN;
        mgm{cell_pair(1),cell_pair(2)} = NaN;
        mps{cell_pair(1),cell_pair(2)} = NaN;
        
    end
    
end

corrs.meanCCGs = mcrosscorrs;
corrs.peakshifts = mps;
corrs.geoMeans = mgm;