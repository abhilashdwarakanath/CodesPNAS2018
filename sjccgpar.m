function [ccSJC] = sjccgpar(params,spikeTimes,cell_pair,jitwin,i,tag)

% This function computes the jitter-corrected cross-correlograms for any
% given cell pair. It needs the cell pair ID, the spike times, a tag for ON
% or OFF phase and the params structure. It depends upon the function
% jitter2.m which has been optimised for use in a parfor loop. 
% The params structure contains the number of samples for the spike trains,
% the maximum lag over which the CCG is computed and various other
% parameters like the number of trials etc. Only the size of the Gaussian
% Smoothing Kernel is hard-coded. Change it according to your needs.
% Abhilash Dwarakanath. Final version - 26.06.2015. MPI for biological
% Cybernetics

for r = 1:size(cell_pair{i},1) % Run through all pairs in a particular distance bin
    
    switch tag
        case 1 %stim on
            shuffle = randperm(params.numtrials);
            fra = zeros(params.fulltime,params.numtrials);
            frb = zeros(params.fulltime,params.numtrials);
            frc = zeros(params.fulltime,params.numtrials);
            frasum = zeros(params.numtrials,1);
            frbsum = zeros(params.numtrials,1);
            frcsum = zeros(params.numtrials,1);
            corr_cc1 = zeros(2*params.jmaxlag+1,params.numtrials);
        case 0 %stim off
            shuffle = randperm(params.numtrials);
            fra = zeros(params.fulltime+30,params.numtrials);
            frb = zeros(params.fulltime+30,params.numtrials);
            frc = zeros(params.fulltime+30,params.numtrials);
            frasum = zeros(params.numtrials,1);
            frbsum = zeros(params.numtrials,1);
            frcsum = zeros(params.numtrials,1);
            corr_cc1 = zeros(2*params.jmaxlag+1,params.numtrials);
    end
    
    switch tag
        
        case 1 %stim on
            
            for j = 1:size(spikeTimes,2)
                st1 = zeros(params.fulltime,1);
                st2 = zeros(params.fulltime,1);
                st3 = zeros(params.fulltime,1);

                if numel(spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}<=params.fulltime))~=0
                    sts1 = (spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}<=params.fulltime));
                    st1(sts1) = 1;
                end
                if numel(spikeTimes{cell_pair{i}(r,2),shuffle(j)}(spikeTimes{cell_pair{i}(r,2),shuffle(j)}<=params.fulltime))~=0
                    sts2 = (spikeTimes{cell_pair{i}(r,2),shuffle(j)}(spikeTimes{cell_pair{i}(r,2),shuffle(j)}<=params.fulltime));
                    st2(sts2) = 1;
                    st3 = jitter2(st2,@randn,jitwin); %
                    st3 = st3';
                end
                spikeTrains{1}(:,j) = st1;
                spikeTrains{2}(:,j) = st2;
                spikeTrains{3}(:,j) = st3; % jittered
               
            end
            
        case 0 %stim off
            
            for j = 1:size(spikeTimes,2)
                st1 = zeros(params.fulltime+30,1);
                st2 = zeros(params.fulltime+30,1);
                st3 = zeros(params.fulltime+30,1);
                if numel(spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}>params.fulltime))~=0
                    sts1 = (spikeTimes{cell_pair{i}(r,1),j}(spikeTimes{cell_pair{i}(r,1),j}>params.fulltime));
                    sts1 = sts1-params.fulltime;
                    st1(sts1) = 1;
                end
                if numel(spikeTimes{cell_pair{i}(r,2),shuffle(j)}(spikeTimes{cell_pair{i}(r,2),shuffle(j)}>params.fulltime))~=0
                    sts2 = (spikeTimes{cell_pair{i}(r,2),shuffle(j)}(spikeTimes{cell_pair{i}(r,2),shuffle(j)}>params.fulltime));
                    sts2 = sts2-params.fulltime;
                    st2(sts2) = 1;
                    st3 = jitter2(st2,@randn,jitwin); % Get Shervin to check this!
                    st3 = st3';
                end
                
                spikeTrains{1}(:,j) = st1;
                spikeTrains{2}(:,j) = st2;
                spikeTrains{3}(:,j) = st3; % jittered

            end
            
    end
     
    % Extract the firing rates of all the trials of the first cell, the
    % second cell and the jittered second cell trials and compute their
    % sum
    
    for tt = 1:params.numtrials
        fra(:,tt) = spikeTrains{1}(:,tt); %Cell 1
        frasum(tt) = sum(fra(:,tt));
        
        frb(:,tt) = spikeTrains{2}(:,tt); % Cell 2
        frbsum(tt) = sum(frb(:,tt));
        
        frc(:,tt) = spikeTrains{3}(:,tt); % Cell 2 jittered
        frcsum(tt) = sum(frc(:,tt));

    end
    
    % Get the mean of their firing rates
    
     persec = params.fulltime/params.Fs;
    
    frassum = (mean(frasum))/persec;
    frbssum = (mean(frbsum))/persec;
    frcssum = (mean(frcsum))/persec;
    tic;
    
    msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair{i},1));
    disp(msg)
    
    % Computer the CCGs
    
    cc1=zeros(2*params.jmaxlag+1,params.numtrials);
    cc2=zeros(2*params.jmaxlag+1,params.numtrials);
    differ=zeros(2*params.jmaxlag+1,params.numtrials);
    
    for tt = 1:params.numtrials
        [cc1(:,tt),lags] = xcorr(fra(:,tt),frb(:,tt),params.jmaxlag); % Expected CCG
        [cc2(:,tt)] = xcorr(fra(:,tt),frc(:,tt),params.jmaxlag); % Jittered CCG
        differ(:,tt) = cc1(:,tt)-cc2(:,tt); % Compute the jitter-corrected CCG i.e. expectedCCG-shuffledCCG
        corr_cc1(:,tt) = differ(:,tt);
    end
    
    el=toc; % elapsed time
    
    msg=sprintf('Time elapsed for pair : %d s',el);
    disp(msg)
    
    % Calculate means of all trials for a given pair
    
    ker = gausswin(5);
    theta = (params.fulltime+1)/params.Fs-abs(lags/params.Fs);
    if isempty(corr_cc1) ~=1
        corr_cc = nanmean2(corr_cc1,2);%sum(corr_cc1,2)./sum(corr_cc1~=0,2);
        corr_cc = conv(corr_cc,ker,'same');
        ccSJC(:,r) = corr_cc/(geomean([frassum frbssum]));
        ccSJC(:,r) = ccSJC(:,r)'./theta;
    else
        ccSJC(:,r) = nan(2*params.jmaxlag+1,1);
    end
    
end

ccSJC = mean(ccSJC,2);
end