function spkTimesSN = cnvrt_spikeMTime2sampleNum(data, favSpikeField, varargin)
% spkTimesSN = cnvrt_spikeMTime2sampleNum(data, favSpikeField, varargin)
% *** NOTE: diff with previous one is the output structure, spkTimesSN{{nUnit}{nTrials} rather spkTimesSN{nUnit, nTrials}  
% e.g. spkTimesSN = cnvrt_SpikeMTimes2sampleNum(data, 'allSpikeTime');
% this function have the spike in times in terms of 'machine times' and
% convert them to corresponding points in term of sample number.
% because after downsampling ...
% It can be useful when you need to find the downsampled index or based
% machene time
% output is a cell array which contain both the spike time (in terms of
% sample number) and also in terms of machene time
% if you look a matrix within one cell (nSpk x 2) the 1st column is in
% sample number the 2nd in machine time.
% *** because 2nd column contain big numbers it's represented in scientific
% notation so the fist column also will be represended like that although
% they are max 4 integer number (don't panic :) )
% ------
% Input:
% 1     data: structure variable which contain data.exp for nChannels and ...
% 2     favSpikeField: indicate the spiking activity (e.g. MUA) in the data 
% (3)   maxNumUnit: the maximum number of channel you want the conversion
%       happen. e.g. you do be done only on first 30 units
%       
% Output:
% 1     spikeTrainS: the spiketrain, which you want to convert its time. It can
%       be part of 'data' (e. g. data.allSpikesByTime) or it can be some external 
%       spikeTrain (e. g. OFFSpikesByTime)
%
% ------
% potential improvments:
% (1) should make some IF weather the input file is cell is mat?
% how many trials how many channels
% 
% ------
% Code Info:
%   creation: some time in 2012 or 2013 by ShS -> shervin.safavi@gmail.com
%   modification: 
%       $ 2014-08-01 by ShS (basic code)
%       $ 2014-08-22 by ShS (compatiblity w/ SUA)
%       $ 2015-05-31 by ShS add an optional input for max number of
%       units/channels
%       $ 2015-06-01 by ShS change the format of output structure

% nTrials = 200; % may be u can find it later in the data file
nTrials = size(data.lfpByTime.t, 2); % number of trials
% nChannels = numel(data.exp.tetids);
% samplingPoints = data.lfp.t{1,1}';
% spikeTrainS = data.allSpikesByTime;
spikeTrainS = data.(favSpikeField);

%% handle optional inputs (varargin):
if nargin > 2 % if the max number of unit was given as an input
    nUnit = varargin{1};
else % if the max number of unit was NOT given as an input
    nUnit = size(spikeTrainS, 1); 
end
% the nUnit can be number of channels, number of single units etc

%%
spkTimesSN{nUnit}{nTrials} = [];
% *** please note that in the original data structure we have
% spkTimes{nUnit,1}{nTrials}; but I don't see the pont for making this
% columnar and the row-wise stucture, so I kept both of them as the default
% orientation of the MATLAB

for iTr = 1 : nTrials
    % explain later about trials
    tmpSamplingPoints = data.lfpByTime.t{1,iTr}';
    for iUnit = 1 : nUnit        
%         tmpRes = data.allSpikesByTime{channelNum,1}{iTr}(:, 1);   
        tmpRes = spikeTrainS{iUnit,1}{iTr}(:, 1); 
        nSpike_tmp = numel(tmpRes);
        DSndx = nan(nSpike_tmp, 1);
        for spkNum = 1 : nSpike_tmp
            tmpDiff = tmpRes(spkNum) - tmpSamplingPoints;
            [~, DSndx(spkNum)] = min(abs(tmpDiff));
        end
        spkTimesSN{iUnit}{iTr}(:, 1) = DSndx;        
    end
end
