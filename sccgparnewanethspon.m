function corrs = sccgparnewanethspon(params,SUspikesByTime,pairs,ts)

tf = 500;
decimFac = tf/params.fs;

es = ts(2)-ts(1);
es = es/1e3;
T = 1/tf:1/tf:es;

ss = ts(1);

% Do first for 90 to 270

for i = 1:size(pairs,1)
    tic;
    cell_pair = pairs(i,:);
    fprintf('Starting a cell pair...\n');
    cell1 = SUspikesByTime{cell_pair(1)};

    cell2 = SUspikesByTime{cell_pair(2)};

            spktimes1 = cell1-ss;
            spktimes1 = ceil((spktimes1.*32.556)*decimFac);
            
            spktimes2 = cell2-ss;
            spktimes2 = ceil((spktimes2.*32.556)*decimFac);
            
            st1 = zeros(1,length(T));
            st2 = zeros(1,length(T));
            st1(spktimes1) = 1;
            st2(spktimes2) = 1;
            
            fr1 = sum(st1)/(length(st1)/tf);
            fr2 = sum(st2)/(length(st2)/tf);
            
            [cc1,lags] = cross_corr(st1,st2,params.maxlag,(1/tf)); % Shuffled CCG
            theta = length(lags)./(length(lags)-abs(lags/tf));
            cc1 = cc1./geomean([fr1 fr2]);
            corr_cc1 = cc1./theta;
            gm =  geomean([fr1 fr2]);
            [~,idx] = max(abs(corr_cc1));
            if idx==1                 timeDiff = [0];             else             timeDiff = lags(idx);             end
            peakShiftstt = timeDiff;

        toc;
        
        mgm{cell_pair(1),cell_pair(2)} = gm;
        mcrosscorrs{cell_pair(1),cell_pair(2)} = corr_cc1;
        mps{cell_pair(1),cell_pair(2)} = peakShiftstt;
    
end

corrs.meanCCGs = mcrosscorrs;
corrs.peakshifts = mps;
corrs.geoMeans = mgm;