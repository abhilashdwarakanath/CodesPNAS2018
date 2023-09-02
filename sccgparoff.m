function [ccSC, rccg, A1, A2] = sccgparoff(spikeTimes,cell_pair,i)

fulltime = 5000; % ON time i.e. the trial time for which each stimulus stays on (in points. 5000 = 10s)
maxlag = 5000; % maximum lag to compute the crosscorrelogram for

for r = 1:size(cell_pair{i},1)
    
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
        fra(:,tt) = spikeTrains{1}([fulltime+1:2*fulltime],tt);
        frasum(tt) = sum(fra(:,tt));
        
        frb(:,tt) = spikeTrains{2}([fulltime+1:2*fulltime],shuffle(tt));
        frbsum(tt) = sum(frb(:,tt));
        
        frc(:,tt) = spikeTrains{2}([fulltime+1:2*fulltime],tt); % Pick out simultaneous spike train from cell 2
        frcsum(tt) = sum(frc(:,tt));
    end
    
    frassum = mean(frasum);
    frbssum = mean(frbsum);
    frcssum = mean(frcsum);
    
    % Get the mean of their firing rates, not the 
    
    tic;
    
    msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair{i},1));
    disp(msg)
    
    % Extract the firing rates of all the trials of the first cell, the
    % second cell and the shuffled second cell trials and compute their
    % sum
    
    for tt = 1:200
        if frasum(tt) > 0 && frbsum(tt) > 0
            aca = xcorr(fra(:,tt),fra(:,tt),maxlag); % ACG of cell 1 (autocorrelogram)
            acb = xcorr(frb(:,tt),frb(:,tt),maxlag); % ACG of cell 2
            [cc1(:,tt)] = xcorr(fra(:,tt),frb(:,tt),maxlag); % Shuffled CCG
            cc1(:,tt) = cc1(:,tt)./([geomean([frassum frbssum])]); % Normalise by geometric mean of firing rates
            [cc2(:,tt)] = xcorr(fra(:,tt),frc(:,tt),maxlag); % Expected CCG
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
    
    %ker = [0.05 0.25 0.40 0.25 0.05]; % From Smith & Summer 2013
    ker = gausswin(15);
    if isempty(corr_cc1) ~=1
        corr_cc = nanmean2(corr_cc1,2);%sum(corr_cc1,2)./sum(corr_cc1~=0,2);
        corr_cc = conv(corr_cc,ker,'same');
        ccSC(:,r) = corr_cc;
    else
        ccSC(:,r) = nan(2*maxlag+1,1);
    end
    
    intwin = [1 5 10 20 50 100 125 250 500]; % intwin*4 = window in ms
    
    if isempty(acorra) ~=1 && isempty(acorrb) ~=1
        acorrA = nanmean2(acorra,2);
        acorrB = nanmean2(acorrb,2);
        mid = fulltime+1;
        for kk = 1:length(intwin)
            a1(kk) = trapz(acorrA(mid-intwin(kk):mid+intwin(kk)));
            a2(kk) = trapz(acorrB(mid-intwin(kk):mid+intwin(kk)));
            rccg(kk,r)= (trapz(ccSC(mid-intwin(kk):mid+intwin(kk),r)))./(sqrt(a1(kk)*a2(kk)));
            A1(kk,r)=a1(kk);A2(kk,r)=a2(kk);
        end
        
    else
        rccg(kk,r) = nan;
        A1(kk,r)=nan;A2(kk,r)=nan;
    end
    
end