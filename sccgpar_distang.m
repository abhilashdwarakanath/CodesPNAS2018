function [rccg] = sccgpar_distang(params,spikeTimes,cell_pair,tag)

% This function computes the shuffle-corrected cross-correlograms for any
% given cell pair. It needs the cell pair ID, the spike times, a tag for ON
% or OFF phase and the params structure.
% The params structure contains the number of samples for the spike trains,
% the maximum lag over which the CCG is computed and various other
% parameters like the number of trials etc. Only the size of the Gaussian
% Smoothing Kernel is hard-coded. Change it according to your needs.
% Abhilash Dwarakanath. Final version - 26.06.2015. MPI for biological
% Cybernetics

for r = 1:size(cell_pair,1)
    
    switch tag
        case 1 %stim on
            shuffle = randperm(params.numtrials);
            fra = zeros(params.fulltime,params.numtrials);
            frb = zeros(params.fulltime,params.numtrials);
            frc = zeros(params.fulltime,params.numtrials);
            acorra = zeros(2*params.maxlag+1,params.numtrials);
            acorrb = zeros(2*params.maxlag+1,params.numtrials);
            frasum = zeros(params.numtrials,1);
            frbsum = zeros(params.numtrials,1);
            frcsum = zeros(params.numtrials,1);
            corr_cc1 = zeros(2*params.maxlag+1,params.numtrials);
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
                if numel(spikeTimes{cell_pair(r,1),j}(spikeTimes{cell_pair(r,1),j}<=params.fulltime))~=0
                    sts1 = (spikeTimes{cell_pair(r,1),j}(spikeTimes{cell_pair(r,1),j}<=params.fulltime));
                    st1(sts1) = 1;
                end
                if numel(spikeTimes{cell_pair(r,2),j}(spikeTimes{cell_pair(r,2),j}<=params.fulltime))~=0
                    sts2 = (spikeTimes{cell_pair(r,2),j}(spikeTimes{cell_pair(r,2),j}<=params.fulltime));
                    st2(sts2) = 1;
                end
                spikeTrains{1}(:,j) = st1;
                spikeTrains{2}(:,j) = st2;
            end
            
        case 0 %stim off
            
            for j = 1:size(spikeTimes,2)
                st1 = zeros(params.fulltime+30,1);
                st2 = zeros(params.fulltime+30,1);
                if numel(spikeTimes{cell_pair(r,1),j}(spikeTimes{cell_pair(r,1),j}>params.fulltime))~=0
                    sts1 = (spikeTimes{cell_pair(r,1),j}(spikeTimes{cell_pair(r,1),j}>params.fulltime));
                    sts1 = sts1-params.fulltime;
                    st1(sts1) = 1;
                end
                if numel(spikeTimes{cell_pair(r,2),j}(spikeTimes{cell_pair(r,2),j}>params.fulltime))~=0
                    sts2 = (spikeTimes{cell_pair(r,2),j}(spikeTimes{cell_pair(r,2),j}>params.fulltime));
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
        fra(:,tt) = spikeTrains{1}(:,tt);
        frasum(tt) = sum(fra(:,tt));
        
        frb(:,tt) = spikeTrains{2}(:,shuffle(tt));
        frbsum(tt) = sum(frb(:,tt));
        
        frc(:,tt) = spikeTrains{2}(:,tt); % Pick out simultaneous spike train from cell 2
        frcsum(tt) = sum(frc(:,tt));
    end
    
    % Get the mean of their firing rates
    
    frassum = mean(frasum);
    frbssum = mean(frbsum);
    frcssum = mean(frcsum);
    
    tic;
    
    msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair,1));
    disp(msg)
    
    % Computer the CCGs
    
    cc1=zeros(2*params.maxlag+1,params.numtrials);
    cc2=zeros(2*params.maxlag+1,params.numtrials);
    differ=zeros(2*params.maxlag+1,params.numtrials);
    
    for tt = 1:params.numtrials
        aca = xcorr(fra(:,tt),fra(:,tt),params.maxlag); % ACG of cell 1 (autocorrelogram)
        acb = xcorr(frb(:,tt),frb(:,tt),params.maxlag); % ACG of cell 2
        [cc1(:,tt)] = xcorr(fra(:,tt),frb(:,tt),params.maxlag); % Shuffled CCG
        cc1(:,tt) = cc1(:,tt)./(geomean([frassum frbssum])); % Normalise by geometric mean of firing rates
        [cc2(:,tt)] = xcorr(fra(:,tt),frc(:,tt),params.maxlag); % Expected CCG
        cc2(:,tt) = cc2(:,tt)./(geomean([frassum frcssum])); % Normalise by geometric mean of firing rates
        differ(:,tt) = cc2(:,tt)-cc1(:,tt); %#ok<*PFOUS> % Compute the shuffle corrected CCG i.e. expectedCCG-shuffledCCG
        corr_cc1(:,tt) = differ(:,tt);
        acorra(:,tt) = aca;
        acorrb(:,tt) = acb;
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
        ccSC(:,r) = nan(2*params.maxlag+1,1);
    end
    
    % params.intwin*4 = window in ms
    
    if isempty(acorra) ~=1 && isempty(acorrb) ~=1
        acorrA = nanmean2(acorra,2);
        acorrB = nanmean2(acorrb,2);
        mid = params.fulltime+1;
            a1 = trapz(acorrA(mid-params.intwin:mid+params.intwin));
            a2 = trapz(acorrB(mid-params.intwin:mid+params.intwin));
            rccg(r)= (trapz(ccSC(mid-params.intwin:mid+params.intwin,r)))./(sqrt(a1*a2));
        
    else
        rccg(r) = nan;
    end
    
end
%addAttachedFiles('local', rccg)
rccg = mean(rccg);