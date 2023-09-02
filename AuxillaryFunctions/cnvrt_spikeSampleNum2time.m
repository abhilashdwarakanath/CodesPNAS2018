function spikesByTime = cnvrt_spikeSampleNum2time(spikesBySampleNum, SF)
% spikeByTime = cnvrt_spikeSampleNum2time(spikeBySampleNum, SF)
% Convert the spikes sample number to spike time based on the sampling
% frequency.
% The output will be based on "mili second". So, the spikesByTime{1}{1}
% indicates that first spike of the first trials happend at
% "spikesByTime{1}{1}"ms
%
% EXAMPLE: -
%
% ------
% Input:
% 1     spikeBySampleNum: cellArray{nUnit}{nTrials}
% 2     obligatoryVariable2: 3*4 array
% 
% Output:
% 1     spikeByTime: cellArray{nUnit}{nTrials}
%       same size as the input
%
% ------
% see also cnvrt_spikeMTime2sampleNum
% ------
% potential improvments:
% (1) check if all the channels does have equal number of trials if that's
% not the case send a warming to the user; check it by unique(nTr)
% ------
% Code Info:
%   creation: 2015-06-01 by ShS (shervin.safavi@gmail.com)
%   modification:
%       $ 201?

%% some parameters which might be needed
nUnit = length(spikesBySampleNum);
singleSampleDuration = (1/SF)*1000; % in mili second

%% conversion
nTr = nan(nUnit,1);
for iUnit = 1 : nUnit
    nTr(iUnit) = length(spikesBySampleNum{iUnit});
    for iTr = 1 : nTr(iUnit)
        spikesByTime{iUnit}{iTr} = spikesBySampleNum{iUnit}{iTr} * singleSampleDuration;
    end    
end