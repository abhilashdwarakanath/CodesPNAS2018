function [results] = sccorrelograms(dist,spikeTimes,cell_pair)


%% Compute

% Computer two pairs of CCGs, one for low firing rate neurons and one for
% high firing rate neurons

fulltime = 5000; % ON time i.e. the trial time for which each stimulus stays on (in points. 5000 = 10s)
maxlag = 5000; % maximum lag to compute the crosscorrelogram for

for i = 1:length(dist) % Run through all distance bins
    
    for r = 1:size(cell_pair{i},1) % Run through all pairs in a particular distance bin
        
        shuffle = randperm(200);
        fra = zeros(fulltime,200);
        frb = zeros(fulltime,200);
        frc = zeros(fulltime,200);
        acorra = zeros(2*maxlag+1,200);
        acorrb = zeros(2*maxlag+1,200);
        frasum = zeros(200,1);
        frbsum = zeros(200,1);
        frcsum = zeros(200,1);
        corr_cc1 = zeros(2*maxlag+1,200);
        
        for j = 1:size(spikeTimes,2)
            st1 = zeros(10300,1);
            st2 = zeros(10300,1);
            if sum(spikeTimes{cell_pair{i}(r,1),j})~=0
                st1(spikeTimes{cell_pair{i}(r,1),j}) = 1;
            end
            if sum(spikeTimes{cell_pair{i}(r,2),j})~=0
                st2(spikeTimes{cell_pair{i}(r,2),j}) = 1;
            end
            spikeTrains{1}(:,j) = st1;
            spikeTrains{2}(:,j) = st2;
        end
        
        for tt = 1:length(shuffle)
            fra(:,tt) = spikeTrains{1}([1:fulltime],tt);
            frasum(tt) = sum(fra(:,tt));
            
            frb(:,tt) = spikeTrains{2}([1:fulltime],shuffle(tt));
            frbsum(tt) = sum(frb(:,tt));
            
            frc(:,tt) = spikeTrains{2}([1:fulltime],tt); % Pick out simultaneous spike train from cell 2
            frcsum(tt) = sum(frc(:,tt));
        end
        
        frassum = sum(frasum);
        frbssum = sum(frbsum);
        frcssum = sum(frcsum);
        
        tic;
        
        msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair{i},1));
        disp(msg)
        
        % Extract the firing rates of all the trials of the first cell, the
        % second cell and the shuffled second cell trials and compute their
        % sum
        
        for tt = 1:200
            if frasum(tt) > 0 && frbsum(tt) > 0 && frcsum(tt) > 0
            aca = xcorr(fra(:,tt),fra(:,tt),maxlag); % ACG of cell 1 (autocorrelogram)
            acb = xcorr(frb(:,tt),frb(:,tt),maxlag); % ACG of cell 2
            [cc1(:,tt)] = xcorr(fra(:,tt),frb(:,tt),maxlag); % Shuffled CCG
            cc1(:,tt) = cc1(:,tt)./([geomean([frassum frbssum])]); % Normalise by geometric mean of firing rates
            [cc2(:,tt), ] = xcorr(fra(:,tt),frc(:,tt),maxlag); % Expected CCG
            cc2(:,tt) = cc2(:,tt)./([geomean([frassum frcssum])]); % Normalise by geometric mean of firing rates
            differ(:,tt) = cc2(:,tt)-cc1(:,tt); %#ok<*PFOUS> % Compute the shuffle corrected CCG i.e. expectedCCG-shuffledCCG
            corr_cc1(:,tt) = [differ(:,tt)];
            acorra(:,tt) = [aca];
            acorrb(:,tt) = [acb];
            else
                corr_cc1(:,tt) = nan(2*maxlag+1,1);
                acorra(:,tt) = nan(2*maxlag+1,1);
                acorrb(:,tt) = nan(2*maxlag+1,1);
            end
        end
        
        el=toc; % elapsed time
        
        msg=sprintf('Time elapsed for pair : %d s',el);
        disp(msg)
        
        % Calculate means of all trials for a given pair
        
        ker = gausswin(5);
        
        if isempty(corr_cc1) ~=1
            corr_cc = nanmean2(corr_cc1,2);%sum(corr_cc1,2)./sum(corr_cc1~=0,2);
            corr_cc = conv(corr_cc,ker,'same');
            ccSC{i}(:,r) = corr_cc;
        else
            ccSC{i}(:,r) = nan(2*maxlag+1,1);
        end
        
        intwin = [1 5 10 20 50 100 500 1000 2500]; % intwin*4 = window in ms
        
        if isempty(acorra) ~=1 && isempty(acorrb) ~=1
            acorrA = nanmean2(acorrb,2);%sum(acorra,2)./sum(acorra~=0,2);
            acorrB = nanmean2(acorrb,2);
            mid = fulltime+1;
            for kk = 1:length(intwin)
                a1(kk) = trapz(acorrA(mid-intwin(kk):mid+intwin(kk)));
                a2(kk) = trapz(acorrB(mid-intwin(kk):mid+intwin(kk)));
                rccg{i}(kk,r)= (trapz(ccSC{i}(mid-intwin(kk):mid+intwin(kk),r)))./(sqrt(a1(kk)*a2(kk)));
            end
            
        else
            rccg{i}(kk,r) = nan;
        end
        
    end
    
end

%% Get results and store in a structure

results.SC = ccSC;
results.Rccg = rccg;

%    results.lags = lags;
results.dist = dist;

end

