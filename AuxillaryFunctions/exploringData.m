%% 01 wanna see there is any underlying stuctre in spike train of different channsl
nTrials = 1;
nChannels = numel(data.exp.tetids);

counter = 0;
for trialNum = 1 : nTrials
    counter = 0;
    for channelNum = 1 : nChannels
        counter = counter + 1;
        % within each trial seperate spike times for ON and OFF stimulus
        allChannelSpikeTimes{channelNum} = data.allSpikesByTime{channelNum,1}{1,trialNum}(:, 2);
        % *this last 1 might be changed depends on which column of time should be choosed for making decision
        % decide any spike in any channel belongs to which category
        
        % find the index of onStimulSpike
        % one more efficient way is find 1st ON spike consider all the
        % spike before that spikes happen in ON part of stimulus and
        % whatever after that OFF. [in this way we'll use less memory]
        
        % put this IF because in the first iteration we dont have any onRes/offRes to
        % be cancatinated on the old one.
        %         if counter == 1
        %             onNdx = find(data.allSpikesByTime{channelNum,1}{1,trialNum}(:, 3) > 0);
        %             onResTrials{trialNum} = tmpRes(onNdx);
        %
        %             offNdx = find(data.allSpikesByTime{channelNum,1}{1,trialNum}(:, 3) < 0);
        %             offResTrials{trialNum} = tmpRes(offNdx);
        %         else
        %             onNdx = find(data.allSpikesByTime{channelNum,1}{1,trialNum}(:, 3) > 0);
        %             onResTrials{trialNum} = [onResTrials{trialNum}; tmpRes(onNdx)];
        %
        %             offNdx = find(data.allSpikesByTime{channelNum,1}{1,trialNum}(:, 3) < 0);
        %             offResTrials{trialNum} = [offResTrials{trialNum}; tmpRes(offNdx)];
        %         end
    end
end


rasterplot(allChannelSpikeTimes)


%% 02 exactly 01, but for all the channels
% wanna see there is any underlying stuctre in spike train of different channels
% so we plot the rasterPlot for all trialls. in each subplot Y axis is
% channels
nTrials = 10;
nChannels = numel(data.exp.tetids);

counter = 0;
for trialNum = 1 : nTrials
    counter = 0;
    clear allChannelSpikeTimes
    for channelNum = 1 : nChannels
        counter = counter + 1;
        allChannelSpikeTimes{channelNum} = data.allSpikesByTime{channelNum,1}{1,trialNum}(:, 1);
    end
    %     subplot(10, 10, trialNum)
    rasterplot(allChannelSpikeTimes)
    axis off; axis tight
    title(num2str(trialNum))
    waitforbuttonpress
end

% for some number of random trials, same in waht we did above
nTrials = 12;
nChannels = numel(data.exp.tetids);

favTrials = randi(100,1, nTrials);
counter = 0;
for trialNum = favTrials
    counter = counter + 1;
    clear allChannelSpikeTimes
    for channelNum = 1 : nChannels
        allChannelSpikeTimes{channelNum} = data.allSpikesByTime{channelNum,1}{1,trialNum}(:, 1);
    end
    subplot(3, 4, counter)
    rasterplot(allChannelSpikeTimes)
    axis off; axis tight
    title(num2str(favTrials(counter)))
end

%% 03 phase of spikes on one channel
% pick one channel and one trials
nTrials = 1;
nChannels = 1;
trialNum = 1;
channelNum = 2;

tmpLFP = data.lfp.v{channelNum, 1, trialNum};
tmpLFPtime = data.lfp.t{1,1}';
tmpRes = data.allSpikesByTime{channelNum,1}{trialNum}(:, 1);

% filter parameters
sampleRate = 32000;
lfpSampleRate = 500;
FS = lfpSampleRate;

% spike times changed because of down sampling
% tmpResDS = round((tmpRes/sampleRate)*lfpSampleRate);

wn = [15 30]/(FS/2);
n = 4; % Order of the filter
[b,a] = butter(n, wn, 'bandpass');

filtLFP_temp = filtfilt(b, a, tmpLFP); %
AnlcLFP_temp = hilbert(filtLFP_temp);
Phase_temp = angle(AnlcLFP_temp);

% let's look the spikes phase histogram
% for doing that first of all we need the spike times in with new sampling
% rate
% nSpike_tmp = numel(tmpRes);
% DSndx = nan(nSpike_tmp, 1);
% for spkNum = 1 : nSpike_tmp
% tmpDiff = tmpRes(spkNum) - tmpLFPtime;
% [~, DSndx(spkNum)] = min(abs(tmpDiff));
% end



% for spkNum = 1 : numel(tmpRes)
%     spkNdx(spkNum) = find(tmpLFPtime == tmpRes(spkNum));
% end
% spkNdx = find(tmpRes == tmpLFPtime)
hist(Phase_temp(allSpikesByIndex{channelNum, trialNum}(:, 1)), 50);


%% 04 very similar to 03 but for random chose of trialNum and channelNum
% it wont conisder trial-channels with less than 100 spikes

for i = 1 : 100
    trialNum = randi(200);
    channelNum = randi(93);
    
    tmpLFP = data.lfp.v{channelNum, 1, trialNum};
    tmpRes = data.allSpikesByTime{channelNum,1}{trialNum}(:, 1);
    % only care when u have many spikes
    if numel(tmpRes) < 800 i = i - 1;  continue; end
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    
    % spike times changed because of down sampling
    % tmpResDS = round((tmpRes/sampleRate)*lfpSampleRate);
    
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    filtLFP_temp = filtfilt(b, a, tmpLFP); %
    AnlcLFP_temp = hilbert(filtLFP_temp);
    Phase_temp = angle(AnlcLFP_temp);
    
    hist(Phase_temp(allSpikesByIndex{channelNum, trialNum}(:, 1)), 80);
    axis tight
    title(['trial ', num2str(trialNum), ' ,channel ', num2str(channelNum)]);
    
    waitforbuttonpress
    
end

%% 05 similar to 04 but cancatinate the spiks from all channels (in one trial) and plot the histogram
for counter = 1 : 20
    % trialNum = 1;
    trialNum = randi(200);
    nChannels = numel(data.exp.tetids);
    allSpikePhase = [];
    for channelNum = 1 : nChannels
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        tmpRes = data.allSpikesByTime{channelNum,1}{trialNum}(:, 1);
        % only care when u have many spikes
        %     if numel(tmpRes) < 800 i = i - 1;  continue; end
        
        % filter parameters
        sampleRate = 32000;
        lfpSampleRate = 500;
        FS = lfpSampleRate;
        
        % spike times changed because of down sampling
        % tmpResDS = round((tmpRes/sampleRate)*lfpSampleRate);
        
        wn = [15 30]/(FS/2);
        n = 4; % Order of the filter
        [b,a] = butter(n, wn, 'bandpass');
        
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        spikePhase_temp = Phase_temp(allSpikesByIndex{channelNum, trialNum}(:, 1));
        allSpikePhase = [allSpikePhase; spikePhase_temp];
    end
    hist(allSpikePhase, 80);
    waitforbuttonpress
end

%% ?? like the resultant length, any prefered phase or not



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ON-OFF Seperated %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% what I did in 04 and 05 was look at spiking-phas histogram of ALL the spikes within one trail.
% I didn't care abput those happen in the ON-stimulus part or OFF-stimulus
% Now I gonna do the same, but seperatly for ON and OFF part
%% 01s very similar to 03 (last part) but for random chose of trialNum and channelNum
% it wont conisder trial-channels with less than 100 spikes

for i = 1 : 1000
    trialNum = randi(200);
    channelNum = randi(93);
    
    tmpLFP = data.lfp.v{channelNum, 1, trialNum};
    tmpRes = data.allSpikesByTime{channelNum,1}{trialNum}(:, 1);
    % only care when u have many spikes
    if numel(tmpRes) < 800 i = i - 1;  continue; end
    disp(i)
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    
    % spike times changed because of down sampling
    % tmpResDS = round((tmpRes/sampleRate)*lfpSampleRate);
    
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    filtLFP_temp = filtfilt(b, a, tmpLFP); %
    AnlcLFP_temp = hilbert(filtLFP_temp);
    Phase_temp = angle(AnlcLFP_temp);
    
    % hist for ON part
    subplot 211
    hist(Phase_temp(ONSpikesByIndex{channelNum, trialNum}(:, 1)), 80);
    title(sprintf(['trial ', num2str(trialNum), ' ,channel ', num2str(channelNum), '\nON']));
    %     axis tight
    
    % hist for OFF part
    subplot 212
    hist(Phase_temp(OFFSpikesByIndex{channelNum, trialNum}(:, 1)), 80);
    title('OFF')
    
    waitforbuttonpress
    
end

%% 02s similar to 01s but cancatinate the spiks from all channels (in one trial) and plot the histogram
for counter = 1 : 20
    % trialNum = 1;
    trialNum = randi(200);
    nChannels = numel(data.exp.tetids);
    %     allSpikePhase = [];
    ONSpikePhase = [];
    OFFSpikePhase = [];
    
    for channelNum = 1 : nChannels
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        tmpRes = data.allSpikesByTime{channelNum,1}{trialNum}(:, 1);
        
        % filter parameters
        sampleRate = 32000;
        lfpSampleRate = 500;
        FS = lfpSampleRate;
        wn = [15 30]/(FS/2);
        n = 4; % Order of the filter
        [b,a] = butter(n, wn, 'bandpass');
        
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        % ON
        spikePhase_temp = Phase_temp(ONSpikesByIndex{channelNum, trialNum}(:, 1));
        ONSpikePhase = [ONSpikePhase; spikePhase_temp];
        
        % OFF
        spikePhase_temp = Phase_temp(OFFSpikesByIndex{channelNum, trialNum}(:, 1));
        OFFSpikePhase = [OFFSpikePhase; spikePhase_temp];
    end
    
    % hist for ON
    subplot 211
    hist(ONSpikePhase, 80);
    title(sprintf(['trial ', num2str(trialNum), '\nON']));
    axis tight
    
    % hist for OFF
    subplot 212
    hist(OFFSpikePhase, 100);
    title('OFF')
    axis tight
    
    waitforbuttonpress
end

%% 03s exactly 02s but different representation
% in this representation we indicate the phase which MIGHT have locking
for counter = 1 : 10
    % trialNum = 1;
    trialNum = randi(200);
    nChannels = numel(data.exp.tetids);
    %     allSpikePhase = [];
    ONSpikePhase = [];
    OFFSpikePhase = [];
    
    for channelNum = 1 : nChannels
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        tmpRes = data.allSpikesByTime{channelNum,1}{trialNum}(:, 1);
        
        % filter parameters
        sampleRate = 32000;
        lfpSampleRate = 500;
        FS = lfpSampleRate;
        wn = [15 30]/(FS/2);
        n = 4; % Order of the filter
        [b,a] = butter(n, wn, 'bandpass');
        
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        % ON
        spikePhase_temp = Phase_temp(ONSpikesByIndex{channelNum, trialNum}(:, 1));
        ONSpikePhase = [ONSpikePhase; spikePhase_temp];
        
        % OFF
        spikePhase_temp = Phase_temp(OFFSpikesByIndex{channelNum, trialNum}(:, 1));
        OFFSpikePhase = [OFFSpikePhase; spikePhase_temp];
    end
    
    % hist for ON
    subplot 211
    [n, xout] = hist(ONSpikePhase, 80);
    [maxF, phMax] = max(n);
    plot(xout*180/pi, n, 'k'); hold on;
    plot([xout(phMax) xout(phMax)]*180/pi, [0 maxF],'LineWidth',3,'Color',[1 0 0]);
    title(sprintf(['trial ', num2str(trialNum), '\nON']));
    axis tight
    
    % hist for OFF
    subplot 212
    [n, xout] = hist(OFFSpikePhase, 100);
    [maxF, phMax] = max(n);
    plot(xout*180/pi, n, 'k'); hold on;
    plot([xout(phMax) xout(phMax)]*180/pi, [0 maxF],'LineWidth',3,'Color',[1 0 0]);
    title('OFF')
    axis tight
    
    waitforbuttonpress
    clf
end


%% 04s resultant length
nTrials = 200;
R = nan(nTrials, 2);
for trialNum = 1 : nTrials
    nChannels = numel(data.exp.tetids);
    %     allSpikePhase = [];
    ONSpikePhase = [];
    OFFSpikePhase = [];
    
    for channelNum = 1 : nChannels
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        tmpRes = data.allSpikesByTime{channelNum,1}{trialNum}(:, 1);
        
        % filter parameters
        sampleRate = 32000;
        lfpSampleRate = 500;
        FS = lfpSampleRate;
        wn = [15 30]/(FS/2);
        n = 4; % Order of the filter
        [b,a] = butter(n, wn, 'bandpass');
        
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        % ON
        spikePhase_temp = Phase_temp(ONSpikesByIndex{channelNum, trialNum}(:, 1));
        ONSpikePhase = [ONSpikePhase; spikePhase_temp];
        
        % OFF
        spikePhase_temp = Phase_temp(OFFSpikesByIndex{channelNum, trialNum}(:, 1));
        OFFSpikePhase = [OFFSpikePhase; spikePhase_temp];
    end
    R(trialNum, 1) = mean(exp(1i*ONSpikePhase)); % ON
    R(trialNum, 2) = mean(exp(1i*OFFSpikePhase)); % ON
end
imagesc(angle(R));
imagesc(abs(R)); colormap gray


%% 05s look at the sitribution of delta-phiS
% ASSUME that in each trial there is one prefered phase (phase locking). If
% this prefred phase is diffferent for ON and OFF stimulus (I called this delta-phi).
% I gonna look at the distribution of this delta-phiS over trials.
% this delta-phi is constant over trials (a distrubution with a peak at some specific delta-phi)
% i. e. always they have say 50deg difference.
% or at each trial they have some random delta-phi (uniformly distributed)
% phase-locked value can be derived from differnet methods (histogram or resaultant length)

deltaPhiS = angle(R(:,1)) - angle(R(:,2));
[n, xout] = hist(deltaPhiS, 70);
[maxF, phMax] = max(n);
plot(xout*180/pi, n, 'k'); hold on;
plot([xout(phMax) xout(phMax)]*180/pi, [0 maxF],'LineWidth',3,'Color',[1 0 0]);

% or this figure
area(xout*180/pi, n, ...
    'FaceColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
    'EdgeColor',[0.800000011920929 0.800000011920929 0.800000011920929]); hold on;
[n, xout] = hist(deltaPhiS, 70);
[maxF, phMax] = max(n);
plot([xout(phMax) xout(phMax)]*180/pi, [0 maxF],'LineWidth',3,'Color',[1 0 0]);
xlabel('\Delta\phi [deg]')
ylabel('Count of Trials')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase Profiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 01p First we should look at pattern of resultant length across the channels
% ALL spikes
% for each trial we have one pattern of resultant length
nChannels = numel(data.exp.tetids);
nTrials = 200;

% ~ memory allocation
oneTrialLFP = nan(nChannels, numel(data.lfp.v{1, 1, 1}));
oneTrialPhase = nan(nChannels, numel(data.lfp.v{1, 1, 1}));
R = nan(nTrials, nChannels);
% oneTriaRes = [];


for trialNum = 1 : nTrials
    oneTriaRes = [];
    % counter = 0;
    % trialNdx = randi(200, 1, 5)
    % for trialNum = trialNdx
    % for i = 1 : 20
    %     trialNum = randi(200);
    % first make a well structured LFP, PhaseLFP, SpikeTimes for that trial
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    for channelNum = 1 : nChannels
        % LFP, PhaseLFP
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        %         oneTrialLFP(channelNum, :) = tmpLFP;
        oneTrialPhase(channelNum, :) = Phase_temp;
        
        % Spike Times
        tmpRes = allSpikesByIndex{channelNum, trialNum}(:, 1);
        oneTriaRes = [oneTriaRes;  tmpRes];
    end
    
    %  now I shoulf calculate the resultant lengths of channels for ONE trials
    
    % one tmpR belongs to resultant length profile of ONE trials
    % i. e. during spiking activity within one trials how each channel modulated
    % abs(tmpR) give the modulation coeffiecient of each channel and
    % angle(tmpR) in which phase they locked (phase locking make sense only when u have high real value (ans(tmpR)))
    tmpR = mean(exp(1i*oneTrialPhase(:, oneTriaRes)), 2);
    R(trialNum, :) = tmpR';
    
    % counter = counter + 1;
    % R(counter, :) = tmpR';
    
    % imagesc(abs(tmpR)); colorbar
    % title(['Trial ', num2str(trialNum), ', nSpikes ', num2str(numel(oneTriaRes))])
    % caxis([0 0.1])
    
    % waitforbuttonpress
    
end
max(max(abs(R)))
figure
imagesc(abs(R)); colorbar
title('ALL Spikes')
ylabel('Trial Number')
xlabel('Channel Number')


%% 02p similar to 01p, but look at ON and OFF part separatly
% First we should look at pattern of resultant length across the channels
% for each trial we have one pattern of resultant length

% ON
nChannels = numel(data.exp.tetids);
nTrials = 200;

% ~ memory allocation
oneTrialLFP = nan(nChannels, numel(data.lfp.v{1, 1, 1}));
oneTrialPhase = nan(nChannels, numel(data.lfp.v{1, 1, 1}));
Ron = nan(nTrials, nChannels);
% oneTriaRes = [];


for trialNum = 1 : nTrials
    oneTriaRes = [];
    
    % first make a well structured LFP, PhaseLFP, SpikeTimes for that trial
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    for channelNum = 1 : nChannels
        % LFP, PhaseLFP
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        %         oneTrialLFP(channelNum, :) = tmpLFP;
        oneTrialPhase(channelNum, :) = Phase_temp;
        
        % Spike Times
        tmpRes = ONSpikesByIndex{channelNum, trialNum}(:, 1);
        oneTriaRes = [oneTriaRes;  tmpRes];
    end
    
    %  now I shoulf calculate the resultant lengths of channels for ONE trials
    
    % one tmpR belongs to resultant length profile of ONE trials
    % i. e. during spiking activity within one trials how each channel modulated
    % abs(tmpR) give the modulation coeffiecient of each channel and
    % angle(tmpR) in which phase they locked (phase locking make sense only when u have high real value (ans(tmpR)))
    tmpR = mean(exp(1i*oneTrialPhase(:, oneTriaRes)), 2);
    Ron(trialNum, :) = tmpR';
    disp(trialNum);
end

max(max(abs(Ron)))
figure
imagesc(abs(Ron)); colorbar
title('ON Spikes')
ylabel('Trial Number')
xlabel('Channel Number')

% -----------------------------------------
% OFF
nChannels = numel(data.exp.tetids);
nTrials = 200;

% ~ memory allocation
oneTrialLFP = nan(nChannels, numel(data.lfp.v{1, 1, 1}));
oneTrialPhase = nan(nChannels, numel(data.lfp.v{1, 1, 1}));
% Ron = nan(nTrials, nChannels);
% oneTriaRes = [];


for trialNum = 1 : nTrials
    oneTriaRes = [];
    
    % first make a well structured LFP, PhaseLFP, SpikeTimes for that trial
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    for channelNum = 1 : nChannels
        % LFP, PhaseLFP
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        %         oneTrialLFP(channelNum, :) = tmpLFP;
        oneTrialPhase(channelNum, :) = Phase_temp;
        
        % Spike Times
        tmpRes = OFFSpikesByIndex{channelNum, trialNum}(:, 1);
        oneTriaRes = [oneTriaRes;  tmpRes];
    end
    
    %  now I shoulf calculate the resultant lengths of channels for ONE trials
    
    % one tmpR belongs to resultant length profile of ONE trials
    % i. e. during spiking activity within one trials how each channel modulated
    % abs(tmpR) give the modulation coeffiecient of each channel and
    % angle(tmpR) in which phase they locked (phase locking make sense only when u have high real value (ans(tmpR)))
    tmpR = mean(exp(1i*oneTrialPhase(:, oneTriaRes)), 2);
    Roff(trialNum, :) = tmpR';
    disp(trialNum);
end

max(max(abs(Roff)))
figure
imagesc(abs(Roff)); colorbar
title('OFF Spikes')
ylabel('Trial Number')
xlabel('Channel Number')


% I gonna see the distribution of resultant lenghtS
figure
subplot 311
hist(abs(R(:)), 100);
subplot 312
hist(abs(Ron(:)), 100);
subplot 313
hist(abs(Roff(:)), 100);

% all figures together
% phases and R
figure

subplot 321
imagesc(abs(R)); colorbar
title('ALL Spikes')

subplot 323
imagesc(abs(Ron)); colorbar
title('ON Spikes')
ylabel('Trial Number')


subplot 325
imagesc(abs(Roff)); colorbar
title('OFF Spikes')
xlabel('Channel Number')

subplot 322
imagesc(angle(R)); colorbar

subplot 324
imagesc(angle(Ron)); colorbar


subplot 326
imagesc(angle(Roff)); colorbar


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Phase Profiles CHANNELWISE %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 01chw it's similar to last part, but, R of each channel calculated by
% considering the spikes only happen in that channel

% ALL ---------------------------------------------

nChannels = numel(data.exp.tetids);
nTrials = 200;

% ~ memory allocation
Rch = nan(nTrials, nChannels);

for trialNum = 1 : nTrials
    
    % first make a well structured LFP, PhaseLFP, SpikeTimes for that trial
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    for channelNum = 1 : nChannels
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        % Spike Times
        tmpRes = OFFSpikesByIndex{channelNum, trialNum}(:, 1);
        
        % calculate the resultant lengths of channels for ONE trials
        tmpR = mean(exp(1i*Phase_temp(tmpRes)));
        Rch(trialNum, channelNum) = tmpR;
    end
    disp(trialNum);
end

max(max(abs(Rch)))
figure
subplot 211
imagesc(abs(Rch)); colorbar
title('ALL Spikes')
ylabel('Trial Number')
subplot 212
imagesc(angle(Rch)); colorbar
xlabel('Channel Number')

% ON ---------------------------------------------
nChannels = numel(data.exp.tetids);
nTrials = 200;

% ~ memory allocation
Rchon = nan(nTrials, nChannels);

for trialNum = 1 : nTrials
    
    % first make a well structured LFP, PhaseLFP, SpikeTimes for that trial
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    for channelNum = 1 : nChannels
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        % Spike Times
        tmpRes = ONSpikesByIndex{channelNum, trialNum}(:, 1);
        
        % calculate the resultant lengths of channels for ONE trials
        tmpR = mean(exp(1i*Phase_temp(tmpRes)));
        Rchon(trialNum, channelNum) = tmpR;
    end
    disp(trialNum);
end

max(max(abs(Rchon)))
figure
subplot 211
imagesc(abs(Rchon)); colorbar
title('ON Spikes')
ylabel('Trial Number')
subplot 212
imagesc(angle(Rchon)); colorbar
xlabel('Channel Number')

% OFF ---------------------------------------------

nChannels = numel(data.exp.tetids);
nTrials = 200;

% ~ memory allocation
Rchoff = nan(nTrials, nChannels);

for trialNum = 1 : nTrials
    
    % first make a well structured LFP, PhaseLFP, SpikeTimes for that trial
    
    % filter parameters
    sampleRate = 32000;
    lfpSampleRate = 500;
    FS = lfpSampleRate;
    wn = [15 30]/(FS/2);
    n = 4; % Order of the filter
    [b,a] = butter(n, wn, 'bandpass');
    
    for channelNum = 1 : nChannels
        tmpLFP = data.lfp.v{channelNum, 1, trialNum};
        filtLFP_temp = filtfilt(b, a, tmpLFP); %
        AnlcLFP_temp = hilbert(filtLFP_temp);
        Phase_temp = angle(AnlcLFP_temp);
        
        % Spike Times
        tmpRes = OFFSpikesByIndex{channelNum, trialNum}(:, 1);
        
        % calculate the resultant lengths of channels for ONE trials
        tmpR = mean(exp(1i*Phase_temp(tmpRes)));
        Rchoff(trialNum, channelNum) = tmpR;
    end
    disp(trialNum);
end

max(max(abs(Rchoff)))
figure
subplot 211
imagesc(abs(Rchoff)); colorbar
title('OFF Spikes')
ylabel('Trial Number')
subplot 212
imagesc(angle(Rchoff)); colorbar
xlabel('Channel Number')


% summerizing

% distribution of resultant lenghtS
figure
subplot 311
hist(abs(Rch(:)), 100);
subplot 312
hist(abs(Rchon(:)), 100);
subplot 313
hist(abs(Rchoff(:)), 100);

% all figures together
% phases and R
figure

subplot 321
imagesc(abs(Rch)); colorbar
title('ALL Spikes')

subplot 323
imagesc(abs(Rchon)); colorbar
title('ON Spikes')
ylabel('Trial Number')


subplot 325
imagesc(abs(Rchoff)); colorbar
title('OFF Spikes')
xlabel('Channel Number')

subplot 322
imagesc(angle(Rch)); colorbar

subplot 324
imagesc(angle(Rchon)); colorbar


subplot 326
imagesc(angle(Rchoff)); colorbar



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spike Triggered Average %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating spike triggered average to chech there is not any problem of
% time alignment

tau = 10;
spkCut = 30;
% spkCut = 15;
nChannels = numel(data.exp.tetids);


% first let's look at the spike trian in the lfp and check how's good the
% chosen tau

clf
subplot 211
chLFP = data.lfp.v{channelNum, 1, trialNum};
tmpRes = allSpikesByIndex{channelNum, trialNum}(:, 1);
plot(tmpLFP); hold on
for i = 1 : numel(tmpRes)
    line([tmpRes(i) tmpRes(i)], [-1.2 -1.1]*1e4,'Color','r')
end
% spkNum = 550;
line([tmpRes(spkNum) tmpRes(spkNum)], [-1.2 -1.1]*1e4,'Color','g', 'LineWidth',3)
plot(tmpRes(spkNum) - tau : tmpRes(spkNum) + tau, tmpstLFPch, 'g', 'LineWidth',2);

subplot 212
tmpstLFPch = chLFP(tmpRes(spkNum) - tau : tmpRes(spkNum) + tau);
hold on
for i = spkNum - 10 : spkNum + 10
    line([tmpRes(i) tmpRes(i)], [-.05 +.05]*1e4,'Color','r')
end
plot(tmpRes(spkNum) - tau : tmpRes(spkNum) + tau, tmpstLFPch, 'LineWidth',2);
line([tmpRes(spkNum) tmpRes(spkNum)], [-.05 +.05]*1e4,'Color','r', 'LineWidth',3)


% calculation for all channels
stAve_LFP = nan(nChannels, 2 * tau + 1);
for channelNum = 1 : nChannels
    % for channelNum = 1 : 2
    %     nTrials = 5;
    nTrials = 200;
    
    sumChstLFPallTrials = zeros(2 * tau + 1, 1);
    for trialNum = 1 : nTrials
        sumChstLFP1Trials = zeros(2 * tau + 1, 1);
        chLFP = data.lfp.v{channelNum, 1, trialNum};
        tmpRes = allSpikesByIndex{channelNum, trialNum}(:, 1);
        if numel(tmpRes) < 2 * spkCut
            nTrials = nTrials - 1;
            continue;
        end
        for spkNum = spkCut : numel(tmpRes) - spkCut
            tmpstLFPch = chLFP(tmpRes(spkNum) - tau : tmpRes(spkNum) + tau);
            sumChstLFP1Trials = sumChstLFP1Trials + tmpstLFPch;
        end
        sumChstLFPallTrials = sumChstLFPallTrials + (1/numel(tmpRes)) * sumChstLFP1Trials;
    end
    stAve_LFP(channelNum, :) = (1 / nTrials) * sumChstLFPallTrials;
end

% plotting
stA_LFPfig = figure('Color',[1 1 1]);
for i = 1 : nChannels
    subplot (10, 10, i)
    plot(stAve_LFP(i,:), 'LineWidth',2,'Color',[0 0 0])
    axis tight
    axis off
    %     ylabel(num2str(aveSpkCount(i)))
    title(num2str(i))
    if isnan(stAve_LFP(i,:))
        text(20, 0.5, num2str(round(aveSpkCount(i))), 'color', 'r', 'horizontalAlignment', 'center');
    else
        text(1,min(stAve_LFP(i,:)),num2str(round(aveSpkCount(i))), ...
            'horizontalAlignment', 'left', 'verticalAlignment', 'top', ...
            'FontWeight', 'bold', 'color', 'b');
    end
end



%%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% =========================================================================
% EVERY THING BEFORE NEED CORRECTION BECASUE OF LFP CHANNELS -> LFP MATRIX
% =========================================================================
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Phase Profiles CHANNELWISE %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% with a solution for low spike count channels %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For calculating R I consider all the spikes in all trials just wight them
% with 1/nSpk_ChannelNum. In this way thos channel which contain a few
% spikes also can be investigated they have locked to some phase or not.
% Becasue small number of spikes in a channel doesn't mean that they are
% not synchronized.

% you dont have to do too many extra calculation just use the already calculated
% Rch, Rchon, Rchoff and do a without-Nan-average across trials
% because some trials have R = NaN, the average will be NaN, although there
% is some trial real R value, but becasue of NaN they wont included. So, I
% should remove them in the sum.

% I add AT to old Rs to say Average across Trials



% RchAT = mean(Rch, 1);
for channelNum = 1 : nChannels
    oneChR = Rch(:, channelNum);
    RVndx = ~isnan(oneChR); % real value index
    XnanOneChR = oneChR(RVndx);
    RchAT(channelNum) = mean(XnanOneChR)
end


% RchonAT = mean(Rchon, 1);
for channelNum = 1 : nChannels
    oneChR = Rchon(:, channelNum);
    RVndx = ~isnan(oneChR); % real value index
    XnanOneChR = oneChR(RVndx);
    RchonAT(channelNum) = mean(XnanOneChR)
end

% RchoffAT = mean(Rchoff, 1)
for channelNum = 1 : nChannels
    oneChR = Rchoff(:, channelNum);
    RVndx = ~isnan(oneChR); % real value index
    XnanOneChR = oneChR(RVndx);
    RchoffAT(channelNum) = mean(XnanOneChR)
end


figure
subplot 511
bar(abs(RchAT), 'FaceColor',[0 0 0])

subplot 512
bar(abs(RchonAT), 'FaceColor',[0 0 0])

subplot 513
bar(abs(RchoffAT), 'FaceColor',[0 0 0])

subplot 514
imagesc(abs(RchAT)); colorbar

subplot 515
imagesc(angle(RchAT)); colorbar

% representation across with trial-wise one
% phases and R
figure

subplot(9,2, [1 3])
imagesc(abs(Rch)); colorbar; box on
title('ALL Spikes')
subplot 925
imagesc(abs(RchAT)); colorbar; axis off; caxis([0 .45])

subplot(9,2, [7 9])
imagesc(abs(Rchon)); colorbar
title('ON Spikes')
ylabel('Trial Number')
subplot (9, 2, 11)
imagesc(abs(RchonAT)); colorbar; axis off; caxis([0 .45])


subplot(9,2, [13 15])
imagesc(abs(Rchoff)); colorbar
title('OFF Spikes')
xlabel('Channel Number')
subplot (9, 2, 17)
imagesc(abs(RchoffAT)); colorbar; axis off; caxis([0 .45])

% phases
subplot(9,2, [2 4])
imagesc(angle(Rch) * 180 / pi); colorbar; box on
subplot 926
imagesc(angle(RchAT) * 180 / pi); colorbar; axis off;

subplot(9,2, [8 10])
imagesc(angle(Rchon) * 180 / pi); colorbar
subplot (9, 2, 12)
imagesc(angle(RchonAT) * 180 / pi); colorbar; axis off;

subplot(9,2, [14 16])
imagesc(angle(Rchoff) * 180 / pi); colorbar
subplot (9, 2, 18)
imagesc(angle(RchoffAT) * 180 / pi); colorbar; axis off;


%% playing around mapped image to the electrod grid
figure
MUAimage(angle(RchAT), channelsMap, MUAmap);


%% some not-systematic exploration
[Rchon RchonAT] = calcR(lfp, ONSpikesByIndex, [15 30], 1, channelsMap, MUAmap);
[Rchoff RchoffAT] = calcR(lfp, OFFSpikesByIndex, [15 30], 1, channelsMap, MUAmap);

diff = Rchon - Rchoff;
figure;
imagesc(abs(diff)); colorbar;

diff = RchonAT - RchoffAT;
figure;
imagesc(abs(diff)); colorbar;
figure
bar(abs(diff))

figure
hist(abs(diff(:)), 100);


%% Exploring in spike counts
% let's look at the chaster plots

%% Look at R at differnet frequencies ranges

% calclation and brief look
% freq band: 15 30
[Rch.f1530.all.trch Rch.f1530.all.AT] = calcR(lfp, allSpikesByIndex, [15 30], 1, channelsMap, MUAmap);
[Rch.f1530.on.trch Rch.f1530.on.AT] = calcR(lfp, ONSpikesByIndex, [15 30], 1, channelsMap, MUAmap);
[Rch.f1530.off.trch Rch.f1530.off.AT] = calcR(lfp, OFFSpikesByIndex, [15 30], 1, channelsMap, MUAmap);

[Rch.f2228.all.trch Rch.f2228.all.AT] = calcR(lfp, allSpikesByIndex, [22 28], 1, channelsMap, MUAmap);
[Rch.f2228.on.trch Rch.f2228.on.AT] = calcR(lfp, ONSpikesByIndex, [22 28], 1, channelsMap, MUAmap);
[Rch.f2228.off.trch Rch.f2228.off.AT] = calcR(lfp, OFFSpikesByIndex, [22 28], 1, channelsMap, MUAmap);

[Rch.f2024.all.trch Rch.f2024.all.AT] = calcR(lfp, allSpikesByIndex, [20 24], 1, channelsMap, MUAmap);
[Rch.f2024.on.trch Rch.f2024.on.AT]  = calcR(lfp, ONSpikesByIndex, [20 24], 1, channelsMap, MUAmap);
[Rch.f2024.off.trch Rch.f2024.off.AT]  = calcR(lfp, OFFSpikesByIndex, [20 24], 1, channelsMap, MUAmap);

% freq band: 15 20
[Rch.f1520.all.trch Rch.f1520.all.AT] = calcR(lfp, allSpikesByIndex, [15 20], 3, [], [], channelsMap, MUAmap, lowFiringChannels);
[Rch.f1520.on.trch Rch.f1520.on.AT] = calcR(lfp, ONSpikesByIndex, [15 20], 1, channelsMap, MUAmap);
[Rch.f1520.off.trch Rch.f1520.off.AT] = calcR(lfp, OFFSpikesByIndex, [15 20], 1, channelsMap, MUAmap);


% plotting
% comparing with AT
calcR(lfp, [], [], 3, Rch.f1530.all.AT, [], channelsMap, MUAmap, lowFiringChannels);
calcR(lfp, [], [], 3, Rch.f1530.on.AT, [], channelsMap, MUAmap, lowFiringChannels);
calcR(lfp, [], [], 3, Rch.f1530.off.AT, [], channelsMap, MUAmap, lowFiringChannels);

calcR(lfp, [], [], 3, Rch.f2228.all.AT, [], channelsMap, MUAmap, lowFiringChannels);
calcR(lfp, [], [], 3, Rch.f2228.on.AT, [], channelsMap, MUAmap, lowFiringChannels);
calcR(lfp, [], [], 3, Rch.f2228.off.AT, [], channelsMap, MUAmap, lowFiringChannels);

calcR(lfp, [], [], 3, Rch.f2024.all.AT, [], channelsMap, MUAmap, lowFiringChannels);
calcR(lfp, [], [], 3, Rch.f2024.on.AT, [], channelsMap, MUAmap, lowFiringChannels);
calcR(lfp, [], [], 3, Rch.f2024.off.AT, [], channelsMap, MUAmap, lowFiringChannels);

% looking at the phase profile in a bunch of trials
calcR(lfp, [], [], 4, [], Rch.f1530.all.trch, channelsMap, MUAmap, lowFiringChannels);

calcR(lfp, [], [], 4, Rch.f1530.on.AT, [], channelsMap, MUAmap, lowFiringChannels);
calcR(lfp, [], [], 4, Rch.f1530.off.AT, [], channelsMap, MUAmap, lowFiringChannels);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Phase Profiles Channelwise %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% After Observing Sort of Travelling Wave in Spikes %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 01trw
% After we observe the phase profile in elecrode array divided to part,
% we should check how much this pattern is significant
% What I'm going to do is, pick up one channel, and look at phase profile
% of the electrode array, and see the phase pattern is consistent accross
% all the spikes.
% filter parameters
freqBand = [15 30];
trialNum = 1;

% ~ memory allocation
[nR, nC, ~] = size(lfp);
phaseAllchOneTrial = nan(nR, nC);

sampleRate = 32000;
lfpSampleRate = 500;
FS = lfpSampleRate;
wn = freqBand/(FS/2);
n = 4; % Order of the filter
[b, a] = butter(n, wn, 'bandpass');

for channelNum = 1 : nChannels
    tmpLFP = lfp(:, channelNum, trialNum);
    filtLFP_temp = filtfilt(b, a, tmpLFP);
    AnlcLFP_temp = hilbert(filtLFP_temp);
    Phase_temp = angle(AnlcLFP_temp);
    phaseAllchOneTrial(:, channelNum) = Phase_temp;
end

% Spike Times
channelNum = 1;
tmpRes = allSpikesByIndex{channelNum, trialNum}(:, 1);

spkCutI = 100;
spkCutF = spkCutI + 101;
aroundEdge = 10;

tmpPartialLFP = lfp(tmpRes(spkCutI) - aroundEdge : tmpRes(spkCutF) + aroundEdge, :, trialNum);
tmpPartialLFPphase = phaseAllchOneTrial(tmpRes(spkCutI) - aroundEdge : tmpRes(spkCutF) + aroundEdge, :);


figure

subplot(12, 10, (101 : 110))
imagescnan(tmpPartialLFPphase');
% caxis([-150 150]);
%
% subplot(12, 10, (111 : 120))
% imagescnan(lfp(tmpRes(spkCutI) - 10 : tmpRes(spkCutF) + 10, :, trialNum)');
% axis off; box off


% for i = spkCutI : spkCutF
% %     line([tmpRes(i) tmpRes(i)], [channelNum - 5 channelNum + 5],'Color','k')
%     line([tmpRes(i) tmpRes(i)], [1 nChannels],'Color','k')
% end
% % indicate which channel u r looking at
% line([0 length(phaseAllchOneTrial(tmpRes(spkCutI) - 10 : tmpRes(spkCutF) + 10, :))], ...
%     [channelNum channelNum],'Color','k')

% spkNum = 550;

hold on
counter = 0;
correctedTmpRes = tmpRes - tmpRes(spkCutI) + aroundEdge;
for spkNum = spkCutI + 1 : spkCutF - 1
    % phase of one channel at the time of a spike happen in one of the channel
    phaseAllChOnespike = phaseAllchOneTrial(tmpRes(spkNum), :) * 180 / pi;
    counter = counter + 1;
    subplot(12, 10, counter)
    MUAimage(phaseAllChOnespike, channelsMap, MUAmap, 1);
    axis off; box off;
    % colorbar;
    % caxis([-150 150]);
    % waitforbuttonpress
    % clf
    
    % line([tmpRes(spkNum) tmpRes(spkNum)], [channelNum - 5 channelNum + 5],'Color','r', 'LineWidth',2)
    plotBackground(tmpPartialLFP, correctedTmpRes, spkCutI, spkCutF, channelNum);
    
    line([correctedTmpRes(spkNum -1) correctedTmpRes(spkNum - 1)], [1 nChannels],'Color','k')
    line([correctedTmpRes(spkNum) correctedTmpRes(spkNum)], [1 nChannels],'Color','k', 'LineWidth', 1.25)
    
    waitforbuttonpress
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spike-Field Coherency %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I'd like to calculate the spike-field coherency with Chronux
% the plotting contains Spk-LFP coherency of ON-OFF-All for each channel (in original grid)

% initialization
% frequncy range
lowFreq = 0; % in Hertz
upFreq = 50; % in Hertz

% frequncy converter
% Becuase my spike times are indicated by sample number not by time in
% second. Therefore I consider the sampling rate (Fs) 1, i.e. I define a
% new unit for time 2msec = 1 new unit of time.
% Frequencies should be converted to new unit of frequencies
lowBand = lowFreq * 2e-3; % lower Band in your new unit
upBand = upFreq * 2e-3; % upper band in your new unit

% prepare the neural data in the format supported by Chronux
% all spikes
spikeTrainS = allSpikesByIndex;
LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;
% ON spikes
spikeTrainS = ONSpikesByIndex;
LFP_chronuxFormat{2} = convert_LFP_chronux(lfp(1:size(lfp, 1)/2, :, :));
ONSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
allTypeSpikeTrains{2} = ONSpikesByIndex_chronuxFormat;
% OFF spikes i'll do it later, it's a bit more complicated than ON because
% the spikes should shift backward

% allocate a color to each condition
% black(ALL), red(ON), blue(OFF)
color = {'k', 'r', 'b'};

% setting function parameters for Chronux
params.Fs = 1; % sampling frequency
params.fpass = [lowBand upBand]; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.er = [1 0.05]; % population error bars



[MUArow MUAcol] = channelLocFinder((1 : nChannels), channelsMap, MUAmap);

figure('Color',[1 1 1]);
for iCh = 1 : nChannels
    %     for iSpikeTrainType = 1 : numel(allTypeSpikeTrains)
    for iSpikeTrainType = 1 : 1
        dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
        dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;
        [C, phi, S12, S1, S2, f_newUnit, zerosp, confC, phistd] = ...
            coherencycpt(dataLfp, dataSpk, params);
        f = f_newUnit * .5e3;
        % title(num2str(iCh))
        %% plotting
        %         if iSpikeTrainType == 1
        subplot2(10, 10, MUArow(iCh), MUAcol(iCh)); % LI = Linear Index
        %         end
        %         hold on
        plotsig(C, confC, 1, f, color{iSpikeTrainType});
        %         plot(C, color{iSpikeTrainType})
        axis off
        axis tight
        % title(num2str(i))
        %         text(1,min(stAve_LFP(iCh,:)),num2str(round(aveSpkCount(iCh))), ...
        %             'horizontalAlignment', 'left', 'verticalAlignment', 'top', ...
        %             'FontWeight', 'bold', 'color', 'b');
        
        
        [MUAnanRow MUAnanCol] = find(isnan(MUAmap));
        for j = 1 : numel(MUAnanRow)
            LI = subplot2(10, 10, MUAnanRow(j), MUAnanCol(j), ...
                'ZColor',[1 1 1],'YColor',[1 1 1], 'XColor',[1 1 1], ...
                'Color',[0 0 0]); % LI = Linear Indexa
            axis tight
        end
        
    end
end

%% a single plot for only one channel for some checking
iCh = 49;
lowFreq = 0; % in Hertz
upFreq = 200; % in Hertz


fpass = [lowFreq upFreq] * 2e-3;
% prepare the neural data in the format supported by Chronux

% all spikes
spikeTrainS = allSpikesByIndex;
LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;

dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;

% setting function parameters for Chronux
params.Fs = 1; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.er = [1 0.05]; % population error bars

[C, phi, S12, S1, S2, f_newUnit, zerosp, confC, phistd] = ...
    coherencycpt(dataLfp, dataSpk, params);
f = f_newUnit * .5e3;

figure
plotsig(C, confC, 1, f, color{iSpikeTrainType});
title(num2str(iCh))

%% convert the spike times to real time, to be sure that there is any mistake any where
iCh = 49;
lowFreq = 0; % in Hertz
upFreq = 200; % in Hertz
fpass = [lowFreq upFreq]


% all spikes
spikeTrainS = allSpikesByIndex;
for iTrial = 1 : nTrials
    spikeTrainS{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
end
LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;

dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.er = [1 0.05]; % population error bars

[C, phi, S12, S1, S2, f, zerosp, confC, phistd] = ...
    coherencycpt(dataLfp, dataSpk, params);
figure
plotsig(C, confC, 1, f, color{iSpikeTrainType});
title(num2str(iCh))

% the results are exactly the same so we can conclude that nothing is wrong
% with your new unit of time

%% Calculating coherency for all channels
% setting parameters
lowFreq = 0; % in Hertz
upFreq = 200; % in Hertz
fpass = [lowFreq upFreq]

% all spikes
spikeTrainS = allSpikesByIndex;

for iCh = 1 : nChannels
    for iTrial = 1 : nTrials
        spikeTrainS{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
    end
    LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
    allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
    allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;
    
    dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
    dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;
    
    % setting function parameters for Chronux
    params.Fs = 500; % sampling frequency
    params.fpass = fpass; % frequency range of interest
    params.tapers = []; % tapers; default -> [3 5]
    params.trialave = 1; % average over trials
    %     params.err = [1 0.05]; % population error bars
    params.err = 0; % population error bars
    
    %     [coh, cohPhi, cohS12, cohS1, cohS2, f, zerosp, cohConfC, cohPhistd] = ...
    %         coherencycpt(dataLfp, dataSpk, params);
    [coh(:, :, iCh), cohPhi(:, :, iCh), cohS12(:, :, iCh), cohS1(:, :, iCh), cohS2(:, :, iCh), cohF, zerosp] = ...
        coherencycpt(dataLfp, dataSpk, params);
end

[MUArow MUAcol] = channelLocFinder((1 : nChannels), channelsMap, MUAmap);

figure('Color',[1 1 1]);
for iCh = 1 : nChannels
    %     for iSpikeTrainType = 1 : numel(allTypeSpikeTrains)
    for iSpikeTrainType = 1 : 1
        %         if iSpikeTrainType == 1
        subplot2(10, 10, MUArow(iCh), MUAcol(iCh)); % LI = Linear Index
        plot(f, 10*log(cohS1(:, :, iCh)), 'k') % LFP
        hold on
        plot(f, 10*log(cohS2(:, :, iCh)), 'r') % Spike
        axis off
        axis tight
    end
end
% for dead channels
[MUAnanRow MUAnanCol] = find(isnan(MUAmap));
for j = 1 : numel(MUAnanRow)
    LI = subplot2(10, 10, MUAnanRow(j), MUAnanCol(j), ...
        'ZColor',[1 1 1],'YColor',[1 1 1], 'XColor',[1 1 1], ...
        'Color',[0 0 0]); % LI = Linear Indexa
    axis tight
end


% let's see them one by one
figure
%%
clf
iCh = 94;
plot(f, 10*log(cohS1(:, :, iCh)), 'k') % LFP
hold on
plot(f, 10*log(cohS2(:, :, iCh)), 'r') % Spike
axis tight
legend('LFP', 'Spike Trian')
xlabel('Frequency [Hz]')
ylabel('Log Amplitude')


% % plotsig(10*log(coh), cohConfC, 1, f);
% % plot(f, coh);
% title(num2str(iCh))
% plot(f, 10*log(cohS2))
% hold on
% plot(f, 10*log(cohS1), 'r')

%% coherogram -> coherency accross time
iCh = 32;
lowFreq = 0; % in Hertz
upFreq = 200; % in Hertz
fpass = [lowFreq upFreq]

movingwin = [200e-3 20e-3];

% convert spike times (wich are based on sample numbers) to second
spikeTrainS = allSpikesByIndex;
for iTrial = 1 : nTrials
    spikeTrainS{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
end
LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;

dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
% params.er = [1 0.05]; % population error bars
params.err = 0; % population error bars

[C,phi,S12,S1,S2,t,f,zerosp] = ...
    cohgramcpt(dataLfp, dataSpk, movingwin,params,0);
figure
imagesc(t,f, C');

% all channels in one plot
lowFreq = 0; % in Hertz
upFreq = 200; % in Hertz
fpass = [lowFreq upFreq];

movingwin = [200e-3 20e-3];

% convert spike times (wich are based on sample numbers) to second
spikeTrainS = allSpikesByIndex;

for iCh = 1 : nChannels
    for iTrial = 1 : nTrials
        spikeTrainS{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
    end
    LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
    allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
    allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;
    
    dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
    dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;
    
    % setting function parameters for Chronux
    params.Fs = 500; % sampling frequency
    params.fpass = fpass; % frequency range of interest
    params.tapers = []; % tapers; default -> [3 5]
    params.trialave = 1; % average over trials
    % params.er = [1 0.05]; % population error bars
    params.err = 0; % population error bars
    
    [C(:, :, iCh), phi(:, :, iCh), S12(:, :, iCh), S1(:, :, iCh), S2(:, :, iCh), t, f, zerosp] = ...
        cohgramcpt(dataLfp, dataSpk, movingwin,params,0);
end

[MUArow MUAcol] = channelLocFinder((1 : nChannels), channelsMap, MUAmap);
figure
for iCh = 1 : nChannels
    subplot2(10, 10, MUArow(iCh), MUAcol(iCh)); % LI = Linear Index
    imagesc(t,f, C(:, :, iCh)');
    axis off
    caxis([0 max(C(:))])
    
    hold on
    line([0 20], [15 15], 'color', 'r')
    line([0 20], [30 30], 'color', 'r')
end

%%
figure

iCh = 66;
imagesc(t,f, C(:, :, iCh)');
caxis([0 max(C(:))])
colorbar
title(num2str(iCh))

hold on
line([0 20], [15 15], 'color', 'r')
line([0 20], [30 30], 'color', 'r')



%% spectral power
% spectral power of spike trian

% spectral power of LFP

% not required because looks decaying very fast

%% spectrogram -> spectrum across time
% It's already calculatedd in calculation of cohergram

% some toy plots for above stuff
[MUArow MUAcol] = channelLocFinder((1 : nChannels), channelsMap, MUAmap);
figure
for iCh = 1 : nChannels
    subplot2(10, 10, MUArow(iCh), MUAcol(iCh)); % LI = Linear Index
    %     imagesc(t,f, S1(:, :, iCh)');
    imagesc(t,f, 10*log(S1(:, :, iCh))');
    axis off
    %     caxis([0 max(C(:))])
    
    hold on
    line([0 20], [15 15], 'color', 'r')
    line([0 20], [30 30], 'color', 'r')
end

iCh = 93;
imagesc(t,f, 10*log(S2(:, :, iCh))');
%     caxis([0 max(C(:))])

hold on
line([0 20], [15 15], 'color', 'r')
line([0 20], [30 30], 'color', 'r')

plot(f, 10*log(S1(1, :, iCh)))


%% spectral and coherency analysis [sumilar to last section]
% in a narrower frequency band
%% Calculating coherency for all channels 10Hz - 55Hz

% setting parameters
lowFreq = 10; % in Hertz
upFreq = 55; % in Hertz
fpass = [lowFreq upFreq]

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars

% all spikes
spikeTrainS = allSpikesByIndex;

for iCh = 1 : nChannels
    %% convert data to Chronux format
    for iTrial = 1 : nTrials
        spikeTrainS{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
    end
    LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
    allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
    allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;
    
    %% temp spike train and LFP
    dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
    dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;
    
    %% calculation
    disp(iCh)
    [coherency.f1055.coh(:, :, iCh), coherency.f1055.cohPhi(:, :, iCh), ...
        coherency.f1055.S12(:, :, iCh), coherency.f1055.S1(:, :, iCh), coherency.f1055.S2(:, :, iCh), ...
        coherency.f1055.f, coherency.f1055.zerosp] = ...
        coherencycpt(dataLfp, dataSpk, params);
end


%% Calculating coherogram for all channels 10Hz - 55Hz

% setting parameters
lowFreq = 10; % in Hertz
upFreq = 55; % in Hertz
fpass = [lowFreq upFreq];

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars

movingwin = [200e-3 20e-3];

% all spikes
spikeTrainS = allSpikesByIndex;


for iCh = 1 : nChannels
    %% convert data to Chronux format
    for iTrial = 1 : nTrials
        spikeTrainS{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
    end
    LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
    allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
    allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;
    
    %% temp spike train and LFP
    dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
    dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;
    
    %% calculation
    disp(iCh)
    [coherogram.f1055.coh(:, :, iCh), coherogram.f1055.cohPhi(:, :, iCh), ...
        coherogram.f1055.S12(:, :, iCh), coherogram.f1055.S1(:, :, iCh), coherogram.f1055.S2(:, :, iCh), ...
        coherogram.f1055.t, coherogram.f1055.f] = ...
        cohgramcpt(dataLfp, dataSpk, movingwin,params,0);
end
%% I ran  similar to former part for coherency amd coherogramfor but 0Hz - 200Hz
% and store them in similar structure
% noticce that I run it from a temporary m file and delete it, so u won't
% find the related code

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0Hz - 200Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Semi Spectral Analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spike-Field Coherency %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%
iCh = 93;

clf

subplot 221 % both LFP & spike spectrume
plot(f, 10*log(cohS1(:, :, iCh)), 'k') % LFP
hold on
plot(f, 10*log(cohS2(:, :, iCh)), 'r') % Spike
axis tight
legend('LFP', 'Spike Trian')
xlabel('Frequency [Hz]')
ylabel('Log Amplitude')
title('LFP & Spike Train Spectrum')

subplot 222 % Spike train spectrogram
imagesc(t,f, 10*log(S2(:, :, iCh))');
%     caxis([0 max(C(:))])
hold on
line([0 20], [15 15], 'color', 'k')
line([0 20], [30 30], 'color', 'k')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title('Spike Train Spectrogram')

subplot 223 % spike-field coherency
plot(cohF, coh(:, :, iCh), 'k')
axis tight
xlabel('Frequency [Hz]')
ylabel('|C|')
title('Spike-Field Coherency')

subplot 224 % spike-filed coherogram
imagesc(t,f, C(:, :, iCh)');
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title('Spike-Filed Coherogram')
hold on
line([0 20], [15 15], 'color', 'k')
line([0 20], [30 30], 'color', 'k')

% CONCLUSION
% favorite coherency peaks in both p12 and p22:
% 11* 17 49 54 55 56 57 58 59 61 62* 63* 66 [70] [71] [72] [73] 74 76 77 78
% 79 80 82 83* 84*

%  need more investigation
% 47 70 71 89 90 93

% peak around 50Hz in LFP spectrume (* -> somehow)
% 1 2 3 4 5 6 7 8 9 10 11 12 13 17 18 19 20 21 22 23 24 25 26 27 28 29 30
% 31 32 33 34 35 36 37 38 39 40 41 42 43 44 47 49* 50 51 52 53 54 55 56 59
% 62* 63* 64* 65 67 68 69 70 71 72 73 74 75 76 77 78* 79 81 82 83 84 85 86
% 87 88 89 90 91* 92* 93

% didn't: 48 57 58 59 50 61 66 70 80

% peak around 150Hz in LFP spectrume
% 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 47 48 49
% 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 71 72 73 74
% 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93

% didn't:

% some specific attributes:

%   2
% p12 some thing specific happen in the time 10s

%   6
% p12 some thing specific happen in the time 10s
% p22 high coherency in favorite band

%   7
% p12 weired discontuinity in spectrugram
% p12 after 10th second some diffenerce in spectrogram
% p21 high coherency in high frequencies

%   8
% p12 weired discontuinity in spectrugram

%   9
% p12 weired discontuinity in spectrugram [may be this becuase
% of the size of sliding window

%   10
% p12 weired discontuinity in spectrugram

%   11
% p12 weired discontuinity in spectrugram
% p12 after 10th second some diffenerce in spectrogram

%   12
% p12 weired discontuinity in spectrugram
% p12 after 10th second some diffenerce in spectrogram

%   13
% p12 weired discontuinity in spectrugram

%   17
% p12 weired discontuinity in spectrugram
% p12 a bit higher frequency in our favorite band

%   19
% p12 spike trian spectrogram
% p12 after 10th second some diffenerce in spectrogram

%   21
% p12 weired discontuinity in spectrugram
% p12 a bit higher frequency in our favorite band

%   23
% p12 spike trian spectrogram
% p12 after 10th second some diffenerce in spectrogram

%   25
% p12 wired spike trian spectrogram

%   30
% p22 a narrow band around 15 Hz in coherogram
% p12 weired discontuinity in spectrugram

%   31
% p12 weired spike train spectrogram
% p12 more power in 2nd 10sec of spike train spectrugram

%   32
% p22 a narrow band around 15 Hz in coherogram
% p12 weired discontuinity in spectrugram
% p12 more power in 2nd 10sec of spike train spectrugram

%   33
% p22 a narrow band around 15 Hz in coherogram
% p12 weired discontuinity in spectrugram
% p12 more power in 2nd 10sec of spike train spectrugram

%   34
% p12 weired discontuinity in spectrugram

%   35
% p12 weired discontuinity in spectrugram
% p22 weired coherogram

%   37
% p12 weired discontuinity in spectrugram

%   38      - low spike count
% p12 weired discontuinity in spectrugram

%   39
% p12 weired discontuinity in spectrugram
% p21 sort of peak in favorite freqiency range

%   40
% p12 spike trian spectrogram

%   41
% p12 spike trian spectrogram

%   42
% p12 weired discontuinity in spectrugram
% p12 more power in 2nd 10sec of spike train spectrugram

%   43
% p12 weired discontuinity in spectrugram
% p21 sort of peak in favorite freqiency range
% p22 a narrow band around 15 Hz in coherogram

%   44
% p12 weired discontuinity in spectrugram
% p22 a narrow band around 15 Hz in coherogram

%   45
% p11 several peaks at 50Hz, 100Hz, 180Hz and 200Hz

%   46
% p11 several peaks at 50Hz, 70Hz, 130Hz and 180Hz

%   48
% p12 weired discontuinity in spectrugram

%   49
% p12 weired discontuinity in spectrugram
% p21 peak around 30Hz
% p22 a distinct band of frequency with high cohernecy

%   51
% p12 weired discontuinity in spectrugram

%   52      - low spike count
% p21 peak around 10Hz

%   54
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 peak around 17Hz
% p22 a distinct band of frequency with high cohernecy

%   55
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p12 more power in 2nd 10sec of spike train spectrugram
% p21 peak around 17Hz
% p22 a distinct band of frequency with high cohernecy

%   56
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 sort of peak around 17Hz
% p22 a distinct band of frequency with high cohernecy

%   57
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 sort of peak around 19Hz-22Hz
% p22 a distinct band of frequency with high cohernecy

%   58
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 sort of peak around 17Hz
% p22 a distinct band of frequency with high cohernecy

%   59
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p21 peak around 17Hz
% p22 a distinct band of frequency with high cohernecy

%   60
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p22 a distinct band of frequency with high cohernecy

%   61
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 peak around 18Hz
% p22 a distinct band of frequency with high cohernecy

%   62
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p21 sort of peak around 19Hz
% p22 a distinct band of frequency with high cohernecy

%   63
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p12 some special after 10sec
% p21 sort of peak around 17Hz
% p22 a distinct band of frequency with high cohernecy.

%   64
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p22 a distinct narrow band of frequency with high cohernecy

%   65
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band

%   66
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p12 some special after 10sec
% p22 sort of distinct narrow band of frequency with high cohernecy.
% p12 after 10th second some diffenerce in spectrogram

%   67
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p22 sort of distinct band of frequency with high cohernecy.
% p12 after 10th second some diffenerce in spectrogram

%   69
% p12 weired discontuinity in spectrugram

%   70
% p12 NOT weired discontuinity in spectrugram
% p12 very  high power in favorite band
% p22 distinct band of frequency with high cohernecy
% p21 very sharp peak around 17Hz

%   71
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p22 VERY distinct band of frequency with high cohernecy
% p21 peak around 15Hz-25Hz

%   72
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p22 VERY distinct band of frequency with high cohernecy
% p21 peak around 16Hz-18Hz

%   73
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p22 VERY distinct band of frequency with high cohernecy
% p21 peak around 16Hz

%   74
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p22 distinct band of frequency with high cohernecy
% p21 peak around 16Hz

%   75
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band

%   76
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band
% p22 sort of distinct band of frequency with high cohernecy
% p21 peak around 17Hz-21Hz

%   77
% p12 a bit weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 peak around 15Hz-20Hz and 26Hz

%   78
% p12 a BIT weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 peak around 18Hz

%   79
% p12 a BIT weired discontuinity in spectrugram
% p12 very high power in favorite band
% p22 distinct band of frequency with high cohernecy (a bit highrt than our favorite band)
% p21 peak around 15Hz-28Hz

%   80
% p12 a BIT weired discontuinity in spectrugram
% p12 very high power in favorite band
% p22 distinct band of frequency with high cohernecy
% p21 peak around 15Hz

%   81
% p12 a BIT weired discontuinity in spectrugram
% p12 very high power in favorite band

%   82
% p12 weired discontuinity in spectrugram
% p12 a bit high power in favorite band
% p22 distinct band of frequency with high cohernecy (a bit highrt than our favorite band)
% p21 peak around 15Hz-28Hz

%   83
% p12 a BIT weired discontuinity in spectrugram
% p12 very  high power in favorite band
% p21 a small peak around 15Hz-28Hz

%   84
% p12 a BIT weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 a small peak around 14Hz-23Hz

%   85
% p12 VERY weired discontuinity in spectrugram
% p21 a small peak around 17Hz-22Hz

%   86
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 a small peak around 17Hz
% p22 sort of distinct band of frequency with high cohernecy (a bit highrt than our favorite band)

% 87
% p12 a BIT weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 sprt of a small peak around 15Hz-35Hz

% 88
% p12 weired discontuinity in spectrugram
% p12 sort of high power in favorite band

% 89
% p12 very weired discontuinity in spectrugram

% 90
% p12 weired discontuinity in spectrugram
% p12 high power in favorite band
% p22 sort of distinct band of frequency with high cohernecy (a bit highrt than our favorite band)

% 92
% p12 a BIT weired discontuinity in spectrugram
% p12 high power in favorite band
% p21 sprt of a small peak around 20Hz-22Hz and 27Hz-28Hz

% 93
% p12 a BIT weired discontuinity in spectrugram
% p12 VERY high power in favorite band
% p21 ???
% p22 sort of distinct band of frequency with LOW  cohernecy

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0Hz - 200Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10Hz - 55Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Semi Spectral Analysis %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spike-Field Coherency %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%
iCh = 1;

clf

subplot 221 % both LFP & spike spectrume
plot(coherency.f1055.f, 10*log(coherency.f1055.S1(:, :, iCh)), 'k') % LFP
hold on
plot(coherency.f1055.f, 10*log(coherency.f1055.S2(:, :, iCh)), 'r') % Spike
axis tight
legend('LFP', 'Spike Trian')
xlabel('Frequency [Hz]')
ylabel('Log Amplitude')
title('LFP & Spike Train Spectrum')

subplot 222 % Spike train spectrogram
imagesc(coherogram.f1055.t, coherogram.f1055.f, 10*log(coherogram.f1055.S2(:, :, iCh))');
%     caxis([0 max(C(:))])
% hold on
% line([0 20], [15 15], 'color', 'k')
% line([0 20], [30 30], 'color', 'k')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title('Spike Train Spectrogram')

subplot 223 % spike-field coherency
plot(coherency.f1055.f, coherency.f1055.coh(:, :, iCh), 'k')
axis tight
xlabel('Frequency [Hz]')
ylabel('|C|')
title('Spike-Field Coherency')

subplot 224 % spike-filed coherogram
imagesc(coherogram.f1055.t, coherogram.f1055.f, coherogram.f1055.coh(:, :, iCh)');
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title('Spike-Filed Coherogram')
% hold on
% line([0 20], [15 15], 'color', 'k')
% line([0 20], [30 30], 'color', 'k')

%% play around ON - OFF spike field coherency
% in a very messy way

%% Calculating coherency for all channels 0Hz - 80Hz

% setting parameters
lowFreq = 0; % in Hertz
upFreq = 80; % in Hertz
fpass = [lowFreq upFreq]

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars

% spikes
% spikeTrainS = allSpikesByIndex;
spikeTrainS = ONSpikesByIndex;

% LFP

lfpOn = LFP(1:5000, :, :);

lfp = lfpOn;


iCh = 57;
%% convert data to Chronux format
for iTrial = 1 : nTrials
    spikeTrainS{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
end
LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;

%% temp spike train and LFP
dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;

%% calculation
[SFNcoherency.f0080.coh, SFNcoherency.f0080.cohPhi, ...
    SFNcoherency.f0080.S12, SFNcoherency.f0080.S1, coherency.f0080.S2, ...
    SFNcoherency.f0080.f, SFNcoherency.f0080.zerosp] = ...
    coherencycpt(dataLfp, dataSpk, params);



plot(coherency.f0080.f, coherency.f0080.coh(:, :, iCh), 'k')
box off
axis tight
xlabel('Frequency [Hz]')
ylabel('Coherence')
title('Spike-Field Coherency')

hold on
tmpH = max(coherency.f0080.coh(:, :, iCh));
hArea = area([10 30], [tmpH tmpH], 'LineStyle','none');

set(hArea, 'FaceColor', [1 0 1])
alpha(.05)

hold on
plot(SFNcoherency.f0080.f, SFNcoherency.f0080.coh, 'r') % ON


% spikes
% spikeTrainS = allSpikesByIndex;
% spikeTrainS = ONSpikesByIndex;
spikeTrainS = OFFSpikesByIndex;

% LFP

lfpOff = LFP(5000:end, :, :);

lfp = lfpOff;


iCh = 57;
%% convert data to Chronux format
for iTrial = 1 : nTrials
    spikeTrainS{iCh, iTrial} = ...
        (allSpikesByIndex{iCh, iTrial}(:, 1) - 5000) * 2e-3;
end
LFP_chronuxFormat{1} = convert_LFP_chronux(lfp);
allSpikesByIndex_chronuxFormat = convert_spike_chronux(spikeTrainS);
allTypeSpikeTrains{1} = allSpikesByIndex_chronuxFormat;

%% temp spike train and LFP
dataLfp = LFP_chronuxFormat{iSpikeTrainType}(iCh).v;
dataSpk = allTypeSpikeTrains{iSpikeTrainType}(iCh).trial;

%% calculation
[SFNcoherency.f0080.coh, SFNcoherency.f0080.cohPhi, ...
    SFNcoherency.f0080.S12, SFNcoherency.f0080.S1, coherency.f0080.S2, ...
    SFNcoherency.f0080.f, SFNcoherency.f0080.zerosp] = ...
    coherencycpt(dataLfp, dataSpk, params);


hold on
plot(SFNcoherency.f0080.f, SFNcoherency.f0080.coh, 'b') % OFF

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Spike-Triggered Covariance %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of LFP Phases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stCov = calcSTcplxCov(phaseLFP, allSpikesByIndex);
% [stCov, STE] = calcSTcplxCov(phaseLFP, allSpikesByIndex);
[~, eigenVectors, eigenValues] = calcSTcplxCov(phaseLFP, allSpikesByIndex, 'phase');

% overlal image from eigen vectors
figure; imagesc(abs(eigenVectors)); axis square
figure; imagesc(angle(eigenVectors)); axis square

subplotInGrid(eigenVectors, 'r', 'Spike-Triggered Covariance of LFP Phase')

% %% substract mean of all channel form each channel ~ wrong algo
% [~, eigenVectors, eigenValues] = calcSTcplxCov(phaseLFP, allSpikesByIndex, 'phase');
% subplotInGrid(eigenVectors, 'r', 'Spike-Triggered Covariance of LFP Phase')
% % u have the save files

%%
iCh = 93;
subplot 121
MUAimage(abs(eigenVectors(iCh, :))); colorbar; axis off; box off;
caxis([min(abs(eigenVectors(:))) max(abs(eigenVectors(:)))])
subplot 122
MUAimage(angle(eigenVectors(iCh, :)) * 180 / pi); colorbar; axis off; box off;
caxis([-150 150])


%% dirty plot
figure
for iCh = 1 : nChannels
    subplot(10,10, iCh)
    imagesc(abs(stCov(:,:,iCh)))
    axis off; box off; %axis square
end

MUAimage(eigenVectors(iCh, :)); axis off; box off;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Spike-Triggered Covariance %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, eigenVectors2, eigenValues2] = calcSTcplxCov(analyticLFP, allSpikesByIndex, 'complex');

% overlal image from eigen vectors
figure; imagesc(abs(eigenVectors2)); axis square
figure; imagesc(angle(eigenVectors2)); axis square

subplotInGrid(eigenVectors2, 'r', 'Spike-Triggered Covariance')



%% to be able to compare eigenVector and eigenVector2 I should normalized
% them:

% compare the absolute values
figure
% subplot 121
imagesc(abs(eigenVectors)/max(abs(eigenVectors(:)))); axis square
colorbar
% subplot 122
figure
imagesc(abs(eigenVectors2)/max(abs(eigenVectors2(:)))); axis square
colorbar
% very differnet amp;itude


% compare the phases values
figure
subplot 121
imagesc(angle(eigenVectors) * 180 / pi); axis square
colorbar
[min(angle(eigenVectors(:))) max(angle(eigenVectors(:)))] * 180 / pi


subplot 122
imagesc(angle(eigenVectors2)* 180 / pi); axis square
colorbar
[min(angle(eigenVectors2(:))) max(angle(eigenVectors2(:)))] * 180 / pi
% the phase pattern is almost the same, but the amplityude not



%% random spike times or something irrelevant to spike times
%    for iCh = 1 : nChannelsCovComp
clc
iCh = 34;
iTrial = 79;
STEphase = [];

% RANDOMSpike Times

% % all spikes of channel-trial
% tmpRes = allSpikesByIndex{iCh, iTrial}(:, 1);

% % part of spike train
% tmpRes = allSpikesByIndex{iCh, iTrial}(30 : end - 30 , 1) + 540;

% % random spike train with random number of spikes
% tmpRes = randi(nLFPsamples, randi(nLFPsamples,1), 1);

% random spike train with predefined (e.g. 100) number of spikes
tmpRes = randi(nLFPsamples, 10, 1);

% % pick portion of LFP (length 300)
% rnd = randi([500 nLFPsamples - 500],1);
% rnd
% tmpRes = rnd : rnd + 300;



numel(tmpRes)
% STE
tmpSTE = phaseLFP(tmpRes, :, iTrial);
STEphase = [STEphase; tmpSTE];
% end
STE = exp(1i * STEphase);
% dimentionality reduction
[~, eigenValuesRND{iCh, 1}, tmpEigenVec, stCovRND(:, :, iCh)] = doPCA(STE, 1);
eigenVectorsRND(iCh, :) = tmpEigenVec';
eigenValuesRND{iCh, 1}(1:3)


% MUAimage(eigenVectorsRND(iCh, :));
MUAimage(angle(eigenVectorsRND(iCh, :)) * 180/pi);
colorbar
caxis([-150 150])



%% stCov & stAve of RANDOM SPIKE TRAINS
% generate a completely random spike trian for all channels and compute the

for iCh = 1 : nChannels
    for iTrial = 1 : nTrials
        rndSpikeTrainS{iCh, iTrial} = randi(nLFPsamples,size(allSpikesByIndex{iCh, iTrial}));
    end
end
[~, eigenVectorsRndST, eigenValuesRndST] = calcSTcplxCov(phaseLFP, rndSpikeTrainS, 'phase');

figure; imagesc(abs(eigenVectors2)); axis square
figure; imagesc(angle(eigenVectors2)); axis square

subplotInGrid(eigenVectorsRndST, 'r', 'RANDOM SPIKE TRAIN -- Spike-Triggered Covariance')

% rasterplot of this random spike trians
iCh = 4;
for iTrial = 1 : nTrials
    tmpRndSpikeTrains{iTrial} = rndSpikeTrainS{iCh, iTrial}(:, 1);
end
rasterplot(tmpRndSpikeTrains, [], 50);
axis tight

% compute stAve for this random spike trian
% spikeTrainS = rndSpikeTrainS;

% ~ memory allocation
rndRm4 = nan(nChannels, nChannels, nTrials);

for iCh = 1 : nChannels
    for iTrial = 1 : nTrials
        tmpRes = rndSpikeTrainS{iCh, iTrial}(:,1);
        Rm4(iCh, :,  iTrial) = mean(exp(1i*phaseLFP(tmpRes, :, iTrial)));
    end
end

rndRm4AT = nanmean(Rm4, 3);

% spiking channels
subplotInGrid(rndRm4AT, 'r', 'RANDOM SPIKE TRAIN -- Spike Triggered Average');

%% eigen analysis of whole signal
% Let's with whole signal instead of some sample points (spike times)
iTrial = 1;
[~, tmpEigenValues, tmpEigenVec] = doPCA(exp(1i * phaseLFP(:,:, iTrial)), 1);
% eigenVectorsRND(iCh, :) = tmpEigenVec';
tmpEigenValues(1:3)

% MUAimage(tmpEigenVec');
MUAimage(angle(tmpEigenVec') * 180/pi);
colorbar

%% rasterplots
% to see how dense is spiking activity, to have some idea to generate
% random spike trians

% rasterplots
iCh = 4;
for iTrial = 1 : nTrials
    %     tmpSpikeTrains{iTrial} = data.allSpikesByTime{iCh}{iTrial}(:, 2);
    tmpSpikeTrains{iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1);
end
rasterplot(tmpSpikeTrains, [], 50);
axis tight


% compare rasterplots of random and real spike trians
figure
% real spike train
subplot 211
rasterplot(tmpSpikeTrains, [], 20); title('real spike train')
axis tight
% random spike trian
subplot 212
rasterplot(tmpRndSpikeTrains, [], 20); title('random spike train')
axis tight

%% chaster plots
clc
iTrial = 6;
clear tmpSpikeTrains
for iCh = 1 : nChannels
    %     tmpSpikeTrains{iTrial} = data.allSpikesByTime{iCh}{iTrial}(:, 2);
    tmpSpikeTrains{iCh} = allSpikesByIndex{iCh, iTrial}(:, 1);
end
rasterplot(tmpSpikeTrains);
axis tight
clc
disp('done')

% 4 trials
clc
for iTrial = 1 : 4;
    for iCh = 1 : nChannels
        %     tmpSpikeTrains{iTrial} = data.allSpikesByTime{iCh}{iTrial}(:, 2);
        tmpSpikeTrains{iCh} = allSpikesByIndex{iCh, iTrial}(:, 1);
    end
    subplot(4, 1, iTrial)
    rasterplot(tmpSpikeTrains);
    axis tight
end
clc
disp('done')

%% Phase Profiles Channelwise
% After Observing Sort of Travelling Wave in Spikes [phase gradiant in phi(R)]

% After we observe the phase profile in elecrode array divided to part,
% we should check how much this pattern is significant
% What I'm going to do is, pick up one channel, and look at phase profile
% of the electrode array, and see the phase pattern is consistent accross
% all the spikes.
% filter parameters

% we will choose a random channel-trial pair of  and look at the phase profiles as
% explained above

% figure
iteration = 1;
favChannels = [57 59 61 64 68 1 36 35 73 76 86 84 2 67 8 7 39 40 77];
for iteration = 1 : 5
    trialNum = randi(nTrials, 1);
    % then will do the hilbert on all the channels in this trial
    
    % pick a random channel, and check there is enough spike in it or not, if
    % not skip this iteration
    %     randomChannel = randi(nChannels, 1);
    rndNdx = randi(numel(favChannels), 1);
    randomChannel = favChannels(rndNdx);
    
    tmpRes = allSpikesByIndex{randomChannel, trialNum}(:, 1);
    
    % check:
    if numel(tmpRes) < 120
        iteration = iteration - 1;
        numel(tmpRes)
        continue;
    end
    
    %     clc
    %     channelNum
    %     trialNum
    %     readableMUAmap
    
    phaseAllchOneTrial = phaseLFP(:, :, trialNum);
    
    % Spike Times
    channelNum = randomChannel;
    
    maxSpkRand = numel(tmpRes) - 101;
    %     randi(maxSpkRand, 1)
    spkCutI = randi(maxSpkRand, 1);
    spkCutF = spkCutI + 101;
    aroundEdge = 10;
    
    tmpPartialLFP = lfp(tmpRes(spkCutI) - aroundEdge : tmpRes(spkCutF) + aroundEdge, :, trialNum);
    tmpPartialLFPphase = phaseAllchOneTrial(tmpRes(spkCutI) - aroundEdge : tmpRes(spkCutF) + aroundEdge, :);
    
    
    clf % clear every subplot
    
    %     subplot(12, 10, (101 : 110))
    %     imagescnan(tmpPartialLFPphase');
    %     title(['Ch ', num2str(channelNum), ', Trial ', num2str(trialNum), ', nSpk ', num2str(numel(tmpRes))])
    %         caxis([-150 150]);
    
    % only for representing on command window
    clc
    iteration
    dispSpkRange = [spkCutI spkCutF]
    dispLFPrange = [tmpRes(spkCutI) - aroundEdge tmpRes(spkCutF) + aroundEdge]
    disp(diff(dispLFPrange) * 2)
    channelNum
    trialNum
    readableMUAmap
    
    hold on
    counter = 0;
    correctedTmpRes = tmpRes - tmpRes(spkCutI) + aroundEdge;
    
    
    for spkNum = spkCutI + 1 : spkCutF - 1
        
        % phase of one channel at the time of a spike happen in one of the channel
        phaseAllChOnespike = phaseAllchOneTrial(tmpRes(spkNum), :) * 180 / pi;
        counter = counter + 1;
        subplot(12, 10, counter)
        MUAimage(phaseAllChOnespike);
        axis off; box off;
        %         colorbar;
        caxis([-150 150]);
        
        % line([tmpRes(spkNum) tmpRes(spkNum)], [channelNum - 5 channelNum + 5],'Color','r', 'LineWidth',2)
        plotBackground(tmpPartialLFP, correctedTmpRes, spkCutI, spkCutF, channelNum, tmpPartialLFPphase, trialNum, tmpRes);
        
        line([correctedTmpRes(spkNum -1) correctedTmpRes(spkNum - 1)], [1 nChannels],'Color','k')
        line([correctedTmpRes(spkNum) correctedTmpRes(spkNum)], [1 nChannels],'Color','k', 'LineWidth', 1.25)
        
        waitforbuttonpress
    end
end


%% random spike trian but without overlapping with realistic spike train
% if I generate a spike trian based on a uniform random number gerenrator
% it will just a nice sample from LFP signal and if sample big enough only
% reflects spatiotemporal pattern what have been observed in LFP already.
% I gonna exclude the real spike times and generate spike train from time
% points which is left.
% It's not the best possible way, but when I don't know the underlying
% structire of spike trian I can not gerneate a spike trian which does not
% have the those structure.

% a toy spike trian
% choose a iCh and iTrial
iCh = 32;
iTrial = 1;
tmpRealRes = allSpikesByIndex{iCh, iTrial}(:, 1);

% % non-overlapping random spike train
% allowedTimes = setdiff((1:nLFPsamples), tmpRes);
% tmpRndResInd = randi(numel(allowedTimes), numel(tmpRealRes), 1);

% tmpRndRes = allowedTimes(tmpRndResInd);

% selective -SUBJECTIVELY- random spike train
minT = 9600;
maxT = 10000;
tmpPartialRes = tmpRealRes(tmpRealRes <= maxT & tmpRealRes >= minT);

tmpRndRes = randi([minT, maxT], numel(tmpPartialRes), 1);


% compute stCov
STEphase = [];

% STE
tmpSTE = phaseLFP(tmpRndRes, :, iTrial);
% tmpSTE = phaseLFP(tmpRealRes, :, iTrial);
STEphase = [STEphase; tmpSTE];
% end
STE = exp(1i * STEphase);
% dimentionality reduction
[~, eigenValuesRND{iCh, 1}, tmpEigenVec, stCovRND(:, :, iCh)] = doPCA(STE, 1);
eigenVectorsRND(iCh, :) = tmpEigenVec';
eigenValuesRND{iCh, 1}(1:3)


% MUAimage(eigenVectorsRND(iCh, :));
clf
subplot 211
MUAimage(angle(eigenVectorsRND(iCh, :)) * 180/pi);
colorbar
% caxis([-150 150])

subplot 212
rasterRes{1} = tmpRealRes;
rasterRes{2} = tmpRndRes;
rasterplot(rasterRes, [],2)
axis tight

numel(tmpRndRes)
numel(tmpRealRes)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Spike-Triggered Averages %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% + SVD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here I won't drop the amplitude to replicate exactly the complex STA
% look at factor analysis


%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Spike-DOUBLE-Triggered Analysis %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new method of isnpired by spike triggered analysis
% sensitive to spike times and LFP phases

% initial cheking -- not the main method

iTrial = 1;
clear tmpSpikeTrains

favCh = readableMUAmap(:, 5);

for i = 1 : numel(favCh)
    %     tmpSpikeTrains{iTrial} = data.allSpikesByTime{iCh}{iTrial}(:, 2);
    tmpSpikeTrains{i} = allSpikesByIndex{favCh(i), iTrial}(:, 1);
end
rasterplot(tmpSpikeTrains);
axis tight

figure('color', 'w')
for iC = 1 : length(readableMUAmap) - 1
    clear tmpSpikeTrains
    favCh = readableMUAmap(:, iC);
    favCh = favCh(~isnan(favCh));
    
    for i = 1 : numel(favCh)
        %     tmpSpikeTrains{iTrial} = data.allSpikesByTime{iCh}{iTrial}(:, 2);
        tmpSpikeTrains{i} = allSpikesByIndex{favCh(i), iTrial}(:, 1);
    end
    
    subplot(3,3, iC)
    rasterplot(tmpSpikeTrains);
    axis tight
    axis off; box off
end


figure('color', 'w')
for iC = 1 : length(readableMUAmap) - 1
    clear tmpSpikeTrains
    favCh = readableMUAmap(iC, :);
    favCh = favCh(~isnan(favCh));
    
    for i = 1 : numel(favCh)
        %     tmpSpikeTrains{iTrial} = data.allSpikesByTime{iCh}{iTrial}(:, 2);
        tmpSpikeTrains{i} = allSpikesByIndex{favCh(i), iTrial}(:, 1);
    end
    
    subplot(3,3, iC)
    rasterplot(tmpSpikeTrains);
    axis tight
    axis off; box off
end


%% ISI distributions bgjgjvkjgjuhvlk CCG
% let's plot the ISI distribution across all the spikes

iTrial = 1;
tmpSpikeTrains{1} = allSpikesByIndex{3, iTrial}(:, 1);
tmpSpikeTrains{2} = allSpikesByIndex{5, iTrial}(:, 1);

CCGstat = CCG(tmpSpikeTrains);

%
% for iB = 1 : length(readableMUAmap)
iB = 3;
clear tmpSpikeTrains
% favCh = readableMUAmap(iB, :);
favCh = readableMUAmap(:, iB);
favCh = favCh(~isnan(favCh));

for i = 1 : numel(favCh)
    tmpSpikeTrains{i} = allSpikesByIndex{favCh(i), iTrial}(:, 1);
end
[CCGstat, edges] = CCG(tmpSpikeTrains);
% end
figure
%
% tmp = 1./(edges*1e-03);
% tmp(~isfinite(tmp)) = 0;
% plot(tmp, squeeze(CCGstat(1, 2, :)), '.');
plot(edges, squeeze(CCGstat(1, 2, :)));
axis tight

% smoothed version
g = gausswin(20); % <-- this value determines the width of the smoothing window
g = g/sum(g);
smoothedPlot = conv(squeeze(CCGstat(10, 1, :)), g, 'same');
clf
% plot(edges / 2, squeeze(CCGstat(1, 2, :)));
% hold on
plot(edges / 2, smoothedPlot, 'g', 'linewidth', 2)
axis tight

%% 1st 8 neighbors
for iCh = 1 : nChannels;
    tmpCCGstat = zeros(9, 9, 161);
    CCGstat{iCh} = zeros(9, 9, 161);
    % iTrial = 25;
    
    % favCh =
    % 1	2 3
    % 4 5 6
    % 7 8 9
    
    [r, c] = channelLocFinder(iCh, channelsMap, MUAmap);
    clear tmpSpikeTrains
    
    % only for sake of cheking to channel configuration
    favChGrid{iCh} = zeros(10,10);
    
    
    for iTrial = 1 : nTrials
        % 1
        if (r-1 > 0 && c-1 > 0 && ~isnan(readableMUAmap(r-1, c-1)))
            favCh(1) = readableMUAmap(r-1, c-1);
            favChGrid{iCh}(r-1, c-1) = favCh(1);
            tmpSpikeTrains{1} = allSpikesByIndex{favCh(1), iTrial}(:, 1);
        else
            favCh(1) = nan;
            %             favChGrid{iCh}(r, c) = favCh(1);
            tmpSpikeTrains{1} = [];
        end
        
        % 2
        if (r-1 > 0 && ~isnan(readableMUAmap(r-1, c)))
            favCh(2) = readableMUAmap(r-1, c);
            favChGrid{iCh}(r-1, c) = favCh(2);
            tmpSpikeTrains{2} = allSpikesByIndex{favCh(2), iTrial}(:, 1);
        else
            favCh(2) = nan;
            %             favChGrid{iCh}(r, c) = favCh(2);
            tmpSpikeTrains{2} = [];
        end
        
        % 3
        if (r-1 > 0 && c+1 <= 10 && ~isnan(readableMUAmap(r-1, c+1)))
            favCh(3) = readableMUAmap(r-1, c+1);
            favChGrid{iCh}(r-1, c+1) =  favCh(3);
            tmpSpikeTrains{3} = allSpikesByIndex{favCh(3), iTrial}(:, 1);
        else
            favCh(3) = nan;
            %             favChGrid{iCh}(r, c) =  favCh(3);
            tmpSpikeTrains{3} = [];
        end
        
        % 4
        if (c-1 > 0 && ~isnan(readableMUAmap(r, c-1)))
            favCh(4) = readableMUAmap(r, c-1);
            favChGrid{iCh}(r, c-1) = favCh(4);
            tmpSpikeTrains{4} = allSpikesByIndex{favCh(4), iTrial}(:, 1);
        else
            favCh(4) = nan;
            %     favChGrid(r, c-1) = favCh(4);
            tmpSpikeTrains{4} = [];
        end
        
        % 5
        favCh(5) = readableMUAmap(r, c);
        favChGrid{iCh}(r, c) = favCh(5);
        tmpSpikeTrains{5} = allSpikesByIndex{favCh(5), iTrial}(:, 1);
        
        % 6
        if (c+1 <= 10 && ~isnan(readableMUAmap(r, c+1)))
            favCh(6) = readableMUAmap(r, c+1);
            favChGrid{iCh}(r, c+1) = favCh(6);
            tmpSpikeTrains{6} = allSpikesByIndex{favCh(6), iTrial}(:, 1);
        else
            favCh(6) = nan;
            %     favChGrid(r, c+1) = favCh(6);
            tmpSpikeTrains{6} = [];
        end
        
        % 7
        if (c-1 > 0 && r+1 <= 10 && ~isnan(readableMUAmap(r+1, c-1)))
            favCh(7) = readableMUAmap(r+1, c-1);
            favChGrid{iCh}(r+1, c-1) = favCh(7);
            tmpSpikeTrains{7} = allSpikesByIndex{favCh(7), iTrial}(:, 1);
        else
            favCh(7) = nan;
            %     favChGrid(r+1, c-1) = favCh(7);
            tmpSpikeTrains{7} = [];
        end
        
        % 8
        if (r+1 <= 10 && ~isnan(readableMUAmap(r+1, c)))
            favCh(8) = readableMUAmap(r+1, c);
            favChGrid{iCh}(r+1, c) = favCh(8);
            tmpSpikeTrains{8} = allSpikesByIndex{favCh(8), iTrial}(:, 1);
        else
            favCh(8) = nan;
            %     favChGrid(r+1, c) = favCh(8);
            tmpSpikeTrains{8} = [];
        end
        
        % 9
        if (r+1 <= 10 && c+1 <= 10 && ~isnan(readableMUAmap(r+1, c+1)))
            favCh(9) = readableMUAmap(r+1, c+1);
            favChGrid{iCh}(r+1, c+1) = favCh(9);
            tmpSpikeTrains{9} = allSpikesByIndex{favCh(9), iTrial}(:, 1);
        else
            favCh(9) = nan;
            %     favChGrid(r+1, c+1) = favCh(9);
            tmpSpikeTrains{9} = [];
        end
        
        [tmpCCGstat, edges] = CCG(tmpSpikeTrains, 'cal');
        CCGstat{iCh} = CCGstat{iCh} + tmpCCGstat;
        
    end
end
%%
    % favCh =
    % 1	2 3
    % 4 5 6
    % 7 8 9
    
iCh = 61;
clc
% readableMUAmap
favChGrid{iCh}

clf
for k = 1 : 9
    subplot(3, 3, k)
    plot(edges*2, squeeze(CCGstat{iCh}(5, k, :)), 'k');
    set(gca,'XDir','reverse');
    hold on
    %     plot([0 0], [0 max(squeeze(CCGstat(5, k, :)))], 'r')
    plot([0 0], ...
        [min(squeeze(CCGstat{iCh}(5, k, :))) max(squeeze(CCGstat{iCh}(5, k, :)))], 'r')
    axis tight
    
    hold on
    tmpH = max(squeeze(CCGstat{iCh}(5, k, :)));
    hArea = area((1./[-30 30])*1000, [tmpH tmpH], 'LineStyle','none');
    
    set(hArea, 'FaceColor', [1 0 1])
    alpha(.05)
    %     axis off
end

% %% a dirty test
% figure
% for k = 1 : 9
%     subplot(3, 3, k)
%     plot(edges*2, squeeze(CCGstat{iCh}(k, 5, :)));
%     hold on
% %     plot([0 0], [0 max(squeeze(CCGstat(k, 5, :)))], 'r')
%     plot([0 0], ...
%         [min(squeeze(CCGstat{iCh}(k, 5, :))) max(squeeze(CCGstat{iCh}(5, k, :)))], 'r')
%          axis tight
%
%     hold on
%     tmpH = max(squeeze(CCGstat{iCh}(k, 5, :)));
%     hArea = area((1./[-30 30])*1000, [tmpH tmpH], 'LineStyle','none');
%
%     set(hArea, 'FaceColor', [1 0 1])
%     alpha(.05)
%     %     axis off
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Coherency Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part of this coherency analysis have been done before
%% field-field coherency wrt a common LFP channel
% for only one pair (one channel and referece)

refCh = 1;
iCh = 58;

lowFreq = 0; % in Hertz
upFreq = 80; % in Hertz
fpass = [lowFreq upFreq];

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars


LFP_chronuxFormat = convert_LFP_chronux(lfp);

% temp spike train and LFP
dataLfp1 = LFP_chronuxFormat(iCh).v;
dataLfp2 = LFP_chronuxFormat(refCh).v;

% calculation
[FFcoherency.f0080.coh, FFcoherency.f0080.cohPhi, ...
    FFcoherency.f0080.S12, FFcoherency.f0080.S1, FFcoherency.f0080.S2, ...
    FFcoherency.f0080.f] = ...
    coherencyc(dataLfp1, dataLfp2, params);

figure
plot(FFcoherency.f0080.f, FFcoherency.f0080.coh, 'r', 'linewidth', 2)
box off
axis tight
xlabel('Frequency [Hz]')
ylabel('Field-Field Coherence ')
% title('Spike-Field Coherency')


%% field-field coherency wrt a common LFP channel
% for all channel not only one (which was former part)

lowFreq = 0; % in Hertz
upFreq = 80; % in Hertz
fpass = [lowFreq upFreq];

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars


LFP_chronuxFormat = convert_LFP_chronux(lfp);


for iCh = 1 : nChannels
    % temp spike train and LFP
    dataLfp1 = LFP_chronuxFormat(iCh).v;
    dataLfp2 = LFP_chronuxFormat(refCh).v;
    
    % calculation
    [coherency(4).coh(:, iCh), coherency(4).cohPhi(:, iCh), ...
        coherency(4).S12(:, iCh), coherency(4).S1(:, iCh), coherency(4).S2(:, iCh), ...
        coherency(4).f] = ...
        coherencyc(dataLfp1, dataLfp2, params);
end

coherency(4).info.freqBand = [0 80];
coherency(4).info.type = 'Field-Field coherency';
coherency(4).info.ref = ['Channel', num2str(refCh)];
coherency(4).info.nChannels = 93;

figure
%

% hold on
figure

subplot 121
plot(coherency(4).f, coherency(4).coh(:, normalChannels))
hold on
plot(coherency(4).f, coherency(4).coh(:, lowFiringChannels), 'k')
box off
axis tight
xlabel('Frequency [Hz]')
ylabel('Field-Field Coherence')

subplot 122
plot(coherency(4).f, coherency(4).coh(:, normalChannels))
box off
axis tight
xlabel('Frequency [Hz]')
ylabel('Field-Field Coherence')

% let's try coherency of crazy channels in comparasion to normal channels
figure
plot(FFcoherency.f0080.f, FFcoherency.f0080.coh(:, normalChannels), 'g'); % normal channels
hold on
plot(FFcoherency.f0080.f, FFcoherency.f0080.coh(:, lowFiringChannels), 'r'); % crazy channels

% let's try LFP of crazy channels in comparasion to normal channels
figure;
clf
iTrial = 4;
plot(lfp(:, normalChannels, iTrial), 'g'); % normal channels
hold on
plot(lfp(:, lowFiringChannels, iTrial), 'r'); % crazy channels
axis tight

% look like that crazy channels recorded with a very low amplitude

%% let's normalized them and try again

maxLFP = repmat(max(abs(lfp), [], 1), nLFPsamples, 1);
% more accurately should be done in some other way, because we have the
% minus amplitude
lfpNorm = lfp ./ maxLFP;
clear maxLFP

figure;
%
clf
iTrial = 70;
plot(lfpNorm(:, normalChannels, iTrial), 'g'); % normal channels
hold on
plot(lfpNorm(:, lowFiringChannels, iTrial), 'r'); % crazy channels
axis tight
%%
% now compute the coherency with normalized LFP
lowFreq = 0; % in Hertz
upFreq = 80; % in Hertz
fpass = [lowFreq upFreq];

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars


LFP_chronuxFormat = convert_LFP_chronux(lfpNorm);


for iCh = 1 : nChannels
    % temp spike train and LFP
    dataLfp1 = LFP_chronuxFormat(iCh).v;
    dataLfp2 = LFP_chronuxFormat(refCh).v;
    
    % calculation
    [FFcoherency.f0080.coh(:, iCh), FFcoherency.f0080.cohPhi(:, iCh), ...
        FFcoherency.f0080.S12(:, iCh), FFcoherency.f0080.S1(:, iCh), FFcoherency.f0080.S2(:, iCh), ...
        FFcoherency.f0080.f] = ...
        coherencyc(dataLfp1, dataLfp2, params);
end

% let's try coherency of crazy channels in comparasion to normal channels
figure
plot(FFcoherency.f0080.f, FFcoherency.f0080.coh(:, normalChannels), 'g'); % normal channels
hold on
plot(FFcoherency.f0080.f, FFcoherency.f0080.coh(:, lowFiringChannels), 'r'); % crazy channels

% %%
% i = 20;
% plot(FFcoherency.f0080.f, FFcoherency.f0080.coh(:, lowFiringChannels(i)), 'r'); % crazy channels
% lowFiringChannels(i)

%% FF: eigen LFP and coherency with other channels
for iTrial = 1 : nTrials;
    [U, D, ~] = svd(lfp(:, : , iTrial), 'econ');
    eigenLFP(:, iTrial) = U(:, 1);
end

% now compute the coherency with normalized LFP
lowFreq = 0; % in Hertz
upFreq = 80; % in Hertz
fpass = [lowFreq upFreq];

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars

% LFP_chronuxFormat = convert_LFP_chronux(lfpNorm);
LFP_chronuxFormat = convert_LFP_chronux(lfp);

dataLfp2 = eigenLFP;
for iCh = 1 : nChannels
    % temp spike train and LFP
    dataLfp1 = LFP_chronuxFormat(iCh).v;
    
    % calculation
    [coherency(3).coh(:, iCh),  coherency(3).cohPhi(:, iCh), ...
         coherency(3).S12(:, iCh), coherency(3).S1(:, iCh),  coherency(3).S2(:, iCh), ...
         coherency(3).f] = ...
        coherencyc(dataLfp1, dataLfp2, params);
end

coherency(3).info.freqBand = [0 80];
coherency(3).info.type = 'Field-Field coherency';
coherency(3).info.ref = 'eigenLFP';
coherency(3).info.nChannels = 93;


% plot eigen LFP
iTrial = 2;
figure;
hold on
plot(lfpNorm(:, normalChannels, iTrial), 'g'); % normal channels
plot(lfpNorm(:, lowFiringChannels, iTrial), 'r'); % crazy channels
% plot( eigenLFP(:, iTrial),'LineWidth',3,'Color',[0 0 0]);
plot( eigenLFP(:, iTrial) / max(eigenLFP(:, iTrial)),'LineWidth',3,'Color',[0 0 0]);
axis tight

% let's try coherency of crazy channels in comparasion to normal channels
figure
plot(coherency(3).f, coherency(3).coh(:, normalChannels), 'g'); % normal channels
hold on
plot(coherency(3).f, coherency(3).coh(:, lowFiringChannels), 'r'); % crazy channels

%
figure
subplot 121
plot(coherency(3).f, coherency(3).coh(:, normalChannels)); % normal channels
hold on
plot(coherency(3).f, coherency(3).coh(:, lowFiringChannels), 'k'); % crazy channels

subplot 122
plot(coherency(3).f, coherency(3).coh(:, normalChannels)); % normal channels
hold on



%% SF: Spike-eigenLFP coherency eigenLFP and spike train of all other channels
for iTrial = 1 : nTrials;
    [U, D, ~] = svd(lfp(:, : , iTrial), 'econ');
    eigenLFP(:, iTrial) = U(:, 1);
end

% now compute the coherency with normalized LFP
lowFreq = 0; % in Hertz
upFreq = 80; % in Hertz
fpass = [lowFreq upFreq];

% setting function parameters for Chronux
params.Fs = 500; % sampling frequency
params.fpass = fpass; % frequency range of interest
params.tapers = []; % tapers; default -> [3 5]
params.trialave = 1; % average over trials
params.err = 0; % population error bars

% convert data to Chronux format
for iCh = 1 : nChannels
    for iTrial = 1 : nTrials
        spikeTrainByTime{iCh, iTrial} = allSpikesByIndex{iCh, iTrial}(:, 1) * 2e-3;
    end
end
spikesByTime_chronuxFormat = convert_spike_chronux(spikeTrainByTime);

dataLfp = eigenLFP;

for iCh = 1 : nChannels;
    dataSpk = spikesByTime_chronuxFormat(iCh).trial;
    
    % computing
    [coherency(2).coh(:, iCh), coherency(2).cohPhi(:, iCh), ...
        coherency(2).S12(:, iCh), coherency(2).S1(:, iCh), coherency(2).S2(:, iCh), ...
        coherency(2).f, coherency(2).zerosp] = ...
        coherencycpt(dataLfp, dataSpk, params);
end

coherency(2).info.freqBand = [0 80];
coherency(2).info.type = 'Spike-Field coherency';
coherency(2).info.ref = 'eigenLFP';
coherency(2).info.nChannels = 93;

% plot(SFcoherency.f0080.f, SFcoherency.f0080.coh(:, iCh), 'linewidth', 2);
figure; 
clf
plot(coherency(2).f, coherency(2).coh(:, normalChannels));
hold on
plot(coherency(2).f, coherency(2).coh(:, lowFiringChannels), 'k');

subplotInGrid(coherency(2), 'linear', 'Spike-Field Coherency wrt Eiegen LFP')


LFP_chronuxFormat = convert_LFP_chronux(lfp);
% what about coherency wrt it own channel
for iCh = 1 : nChannels;
    dataSpk = spikesByTime_chronuxFormat(iCh).trial;
    dataLfp = LFP_chronuxFormat(iCh).v;
    
    % computing
    [coherency(1).coh(:, iCh), coherency(2).cohPhi(:, iCh), ...
        coherency(2).S12(:, iCh), coherency(2).S1(:, iCh), coherency(2).S2(:, iCh), ...
        coherency(2).f, coherency(2).zerosp] = ...
        coherencycpt(dataLfp, dataSpk, params);
end

coherency(1).info.freqBand = [0 80];
coherency(1).info.type = 'Spike-Field coherency';
coherency(1).info.ref = 'local';
coherency(1).info.nChannels = 93;

figure; 
clf
plot(coherency(1).f, coherency(1).coh(:, normalChannels));
hold on
plot(coherency(1).f, coherency(1).coh(:, lowFiringChannels), 'k');

subplotInGrid(coherency(2), 'linear', 'Spike-Field Coherency wrt Eiegen LFP')
