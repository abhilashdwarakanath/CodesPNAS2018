function [raw,ccSC, rccg, rccg_trialbins] = sccgpar(params,spikeTimes,cell_pair,i,tag)

% This function computes the shuffle-corrected cross-correlograms for any
% given cell pair. It needs the cell pair ID, the spike times, a tag for ON
% or OFF phase and the params structure.
% The params structure contains the number of samples for the spike trains,
% the maximum lag over which the CCG is computed and various other
% parameters like the number of trials etc. Only the size of the Gaussian
% Smoothing Kernel is hard-coded. Change it according to your needs.
% Abhilash Dwarakanath. 

for r = 1:size(cell_pair{i},1)
    
    switch tag
        case 1 %stim on
            shuffle = randperm(params.numtrials); %Shuffle trial order to compute shuffled CCGs
			
			% Initialise all required vectors
			
            fra = zeros(params.fulltime,params.numtrials); % PSTH cell 1	
            frb = zeros(params.fulltime,params.numtrials); % PSTH cell2 SHUFFLED
            frc = zeros(params.fulltime,params.numtrials); % PSTH Cell 2 unshuffled
            acorra = zeros(2*params.maxlag+1,params.numtrials); % Autocorrelogram Cell 1
            acorrb = zeros(2*params.maxlag+1,params.numtrials); % Autocorrelogram Cell 2
            frasum = zeros(params.numtrials,1); % Firing rate
            frbsum = zeros(params.numtrials,1);
            frcsum = zeros(params.numtrials,1);
            corr_cc1 = zeros(2*params.maxlag+1,params.numtrials); % Corrected CCG
        case 0 % stim off
            shuffle = randperm(params.numtrials);
            fra = zeros(params.fulltime+30,params.numtrials);
            frb = zeros(params.fulltime+30,params.numtrials);
            frc = zeros(params.fulltime+30,params.numtrials);
            acorra = zeros(2*params.maxlag+1,params.numtrials);
            acorrb = zeros(2*params.maxlag+1,params.numtrials);
            frasum = zeros(params.numtrials,1);
            frbsum = zeros(params.numtrials,1);
            frcsum = zeros(params.numtrials,1);
            corr_cc1 = zeros(2*params.maxlag+1,params.numtrials);
    end
    
    switch tag
        
        case 1 %stim on
            
            for j = 1:size(spikeTimes,2)
                st1 = zeros(params.fulltime,1);
                st2 = zeros(params.fulltime,1);
                if numel(spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}<=params.fulltime))~=0
                    sts1 = (spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}<=params.fulltime));
                    st1(sts1) = 1; % Create spike trains
                end
                if numel(spikeTimes{cell_pair{i}(r,2),j}(spikeTimes{cell_pair{i}(r,2),j}<=params.fulltime))~=0
                    sts2 = (spikeTimes{cell_pair{i}(r,2),j}(spikeTimes{cell_pair{i}(r,2),j}<=params.fulltime));
                    st2(sts2) = 1;
                end
                spikeTrains{1}(:,j) = st1;
                spikeTrains{2}(:,j) = st2;
            end
            
        case 0 %stim off
            
            for j = 1:size(spikeTimes,2)
                st1 = zeros(params.fulltime+30,1);
                st2 = zeros(params.fulltime+30,1);
                if numel(spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}>params.fulltime))~=0
                    sts1 = (spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}>params.fulltime));
                    sts1 = sts1-params.fulltime;
                    st1(sts1) = 1;
                end
                if numel(spikeTimes{cell_pair{i}(r,2),j}(spikeTimes{cell_pair{i}(r,2),j}>params.fulltime))~=0
                    sts2 = (spikeTimes{cell_pair{i}(r,2),j}(spikeTimes{cell_pair{i}(r,2),j}>params.fulltime));
                    sts2 = sts2-params.fulltime;
                    st2(sts2) = 1;
                end
                spikeTrains{1}(:,j) = st1;
                spikeTrains{2}(:,j) = st2;
            end
            
    end
    
    % Extract the firing rates of all the trials of the first cell, the
    % second cell and the shuffled second cell trials and compute their
    % sum
    
    for tt = 1:length(shuffle)
        fra(:,tt) = (spikeTrains{1}(:,tt));
        frasum(tt) = sum(fra(:,tt));
        
        frb(:,tt) = (spikeTrains{2}(:,shuffle(tt)));
        frbsum(tt) = sum(frb(:,tt));
        
        frc(:,tt) = (spikeTrains{2}(:,tt)); % Pick out simultaneous spike train from cell 2
        frcsum(tt) = sum(frc(:,tt));
    end
    
    % Get the mean of their firing rates
    
    persec = params.fulltime/params.Fs;
    
    frassum = (mean(frasum))/persec;
    frbssum = (mean(frbsum))/persec;
    
    tic;
    
    msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair{i},1));
    disp(msg)
    
    % Computer the CCGs
    
    cc1=zeros(2*params.maxlag+1,params.numtrials);
    cc2=zeros(2*params.maxlag+1,params.numtrials);
    differ=zeros(2*params.maxlag+1,params.numtrials);
    
    %[~,ker1] = smoothingkernel(params.fulltime/params.Fs,params.Fs,0.01,'gaussian'); % 100 ms kernel
    
    for tt = 1:params.numtrials
        %x = conv(fra(:,tt),ker1,'same');
        %y= conv(frb(:,tt),ker1,'same');
        %z = conv(frc(:,tt),ker1,'same');
        %         aca = (xcorr(fra(:,tt),fra(:,tt),params.maxlag)); % ACG of cell 1 (autocorrelogram)
        %         acb = (xcorr(frc(:,tt),frc(:,tt),params.maxlag)); % ACG of cell 2
        %         [cc1(:,tt),lags] = xcorr(fra(:,tt),frb(:,tt),params.maxlag); % Shuffled CCG
        %         [cc2(:,tt)] = xcorr(fra(:,tt),frc(:,tt),params.maxlag); % Expected CCG
        %         differ(:,tt) = cc2(:,tt)-cc1(:,tt); %#ok<*PFOUS> % Compute the shuffle corrected CCG i.e. rawCCG-shuffledCCG
        
        aca = cross_corr(fra(:,tt),fra(:,tt),params.maxlag,1); % ACG of cell 1 (autocorrelogram)
        acb = cross_corr(frc(:,tt),frc(:,tt),params.maxlag,1); % ACG of cell 2
        [cc1(:,tt),lags] = cross_corr(fra(:,tt),frb(:,tt),params.maxlag,1); % Shuffled CCG
        [cc2(:,tt)] = cross_corr(fra(:,tt),frc(:,tt),params.maxlag,1); % Expected CCG
        differ(:,tt) = cc2(:,tt)-cc1(:,tt); %#ok<*PFOUS>
        corr_cc1(:,tt) = differ(:,tt);
        acorra(:,tt) = aca;
        acorrb(:,tt) = acb;
    end
    
    el=toc; % elapsed time
    
    msg=sprintf('Time elapsed for pair : %d s',el);
    disp(msg)
    
    % Collect all trials for a pair
    
    pwsccgs(r,:,:) = corr_cc1;
    accgs1(r,:,:) = acorra;
    accgs2(r,:,:) = acorrb;
    
    % Calculate means of all trials for a given pair
    
    theta = (params.fulltime+1)/params.Fs-abs(lags/params.Fs);
    
    ker = gausswin(5);
    if isempty(corr_cc1) ~=1
        corr_cc = nanmean2(corr_cc1,2);%sum(corr_cc1,2)./sum(corr_cc1~=0,2);
        %corr_cc = conv(corr_cc,ker,'same');
        ccSC(:,r) = corr_cc/(geomean([frassum frbssum])); % Normalise by the geometric mean
        ccSC(:,r) = ccSC(:,r)'./theta; %Normalise by the theta function
        %ccSC(:,r) = corr_cc'./(frassum*frbssum*theta);
        acg1(:,r) = nanmean2(acorra,2)/frassum;
        acg2(:,r) = nanmean2(acorrb,2)/frbssum;
    else
        ccSC(:,r) = nan(2*params.maxlag+1,1);
        acg1(:,r) = nan(2*params.maxlag+1,1);
        acg2(:,r) = nan(2*params.maxlag+1,1);
    end
    
    % Compute mean CCGs for trial bins
    numbins = 11;
    binedges = 1:params.numtrials/(numbins-1):params.numtrials+1;
    
    mid = (params.maxlag)+1;
    rccgtb = zeros(1,numbins-1);
    for n = 1:numbins-1
        ccg_tb = (mean(corr_cc1(:,binedges(n):binedges(n+1)-1),2))/(geomean([frassum frbssum]));
        acorratb = mean(acorra(:,binedges(n):binedges(n+1)-1),2)/frassum;
        a1 = trapz(acorratb);
        acorrbtb = mean(acorrb(:,binedges(n):binedges(n+1)-1),2)/frbssum;
        a2 = trapz(acorrbtb);
        rccgtb(n)= (trapz(ccg_tb(mid-params.intwin(end):mid+params.intwin(end))))/(sqrt(a1*a2));
    end
    
    rccg_trialbins(:,r) = rccgtb;
    
    clear theta
    
    % params.intwin*4 = window in ms
    
    for kk = 1:length(params.intwin)
        %a1(kk) = trapz(acg1(:,r)); % changed to Alex Ecker's method
        %a2(kk) = trapz(acg2(:,r));
        a1(kk) = trapz(acg1(mid-params.intwin(kk):mid+params.intwin(kk),r));
        a2(kk) = trapz(acg2(mid-params.intwin(kk):mid+params.intwin(kk),r));
        rccg(kk,r)= (trapz(ccSC(mid-params.intwin(kk):mid+params.intwin(kk),r)))./(sqrt(a1(kk)*a2(kk)));
        A1(kk,r)=a1(kk);A2(kk,r)=a2(kk);
    end
    
end

ccSC = mean(ccSC,2);
rccg = mean(rccg,2);
rccg_trialbins = mean(rccg_trialbins,2);
raw.sccg = pwsccgs;
raw.accg.a = accgs1;
raw.accg.b = accgs2;