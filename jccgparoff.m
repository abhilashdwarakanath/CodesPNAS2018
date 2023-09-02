function [ccJC] = jccgparoff(spikeTimes,cell_pair,jitwin,i)

fulltime = 5000; % ON time i.e. the trial time for which each stimulus stays on (in points. 5000 = 10s)
maxlag = 1000; % maximum lag to compute the crosscorrelogram for

for r = 1:size(cell_pair{i},1) % Run through all pairs in a particular distance bin
    
    shuffle = randperm(200);
    fra = zeros(fulltime,200);
    frb = zeros(fulltime,200);
    frc = zeros(fulltime,200);
    frasum = zeros(200,1);
    frbsum = zeros(200,1);
    frcsum = zeros(200,1);
    corr_cc1 = zeros(2*maxlag+1,200);
    
    for j = 1:size(spikeTimes,2)
        st1 = zeros(10300,1);
        st2 = zeros(10300,1);
        st3 = zeros(10300,1);
        
        if sum(spikeTimes{cell_pair{i}(r,1),j})~=0
            st1(spikeTimes{cell_pair{i}(r,1),j}) = 1;
        end
        if sum(spikeTimes{cell_pair{i}(r,2),j})~=0
            st2(spikeTimes{cell_pair{i}(r,2),j}) = 1;
            st3 = jitter2(st2,@randn,jitwin); % Get Shervin to check this!
            st3 = st3';
        end
        
        spikeTrains{1}(:,j) = st1;
        spikeTrains{2}(:,j) = st2;
        spikeTrains{3}(:,j) = st3;
    end
    
    for tt = 1:length(shuffle)
        fra(:,tt) = spikeTrains{1}([fulltime+1:2*fulltime],tt);
        frasum(tt) = sum(fra(:,tt));
        
        frb(:,tt) = spikeTrains{2}([fulltime+1:2*fulltime],tt);
        frbsum(tt) = sum(frb(:,tt));
        
        frc(:,tt) = spikeTrains{3}([fulltime+1:2*fulltime],tt);
        frcsum(tt) = sum(frc(:,tt));
        
    end
    
    frassum = mean(frasum);
    frbssum = mean(frbsum);
    frcssum = mean(frcsum);
    
    tic;
    
    msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair{i},1));
    disp(msg)
    
    % Extract the firing rates of all the trials of the first cell, the
    % second cell and the jittered second cell trials and compute their
    % sum
    
    for tt = 1:200
        if frasum(tt) > 0 && frbsum(tt) > 0 && frcsum(tt) > 0
            [cc1(:,tt)] = xcorr(fra(:,tt),frb(:,tt),maxlag); % Expected CCG
            cc1(:,tt) = cc1(:,tt)./([geomean([frassum frbssum])]); % Normalise by geometric mean of firing rates
            [cc2(:,tt)] = xcorr(fra(:,tt),frc(:,tt),maxlag); % Jittered CCG
            cc2(:,tt) = cc2(:,tt)./([geomean([frassum frcssum])]); % Normalise by geometric mean of firing rates
            differ(:,tt) = cc2(:,tt)-cc1(:,tt); % Compute the jitter-corrected CCG i.e. expectedCCG-shuffledCCG
            corr_cc1(:,tt) = [differ(:,tt)];
        else
            corr_cc1(:,tt) = nan(2*maxlag+1,1);
        end
    end
    
    el=toc; % elapsed time
    
    msg=sprintf('Time elapsed for pair : %d s',el);
    disp(msg)
    
    % Calculate means of all trials for a given pair
    
    ker = gausswin(20);
    
    if isempty(corr_cc1) ~=1
        corr_cc = nanmean2(corr_cc1,2);%sum(corr_cc1,2)./sum(corr_cc1~=0,2);
        corr_cc = conv(corr_cc,ker,'same');
        ccJC(:,r) = corr_cc;
    else
        ccJC(:,r) = nan(2*maxlag+1,1);
    end
    
end