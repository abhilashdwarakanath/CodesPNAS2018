function [spikeCounts, nBin] = estimateFiringRate(spikeTrains, binSize, times)
% [spikeCounts, nBin] = estimateFiringRate(spikeTrains, binSize, times)
% This function will find the spike counts within each time-bin
% You can either deal with time or samples you just need to be consistent
% your input
% 
% ------
% Input:
% 1      spikeTrains: a cell array with dim of {nChannel x nTrial} 
% 2      binSize: the favorite size of the bin based on sample size unit 
% 3      times: a structure which should contain the period in which you
%        like to estimate the firing rates. The structure shoud have this
%        fileds: "times.startT" and "times.stopT"
% 
% Output:
% 1      spikeCount: spikes in selected bin size in matrix of nBin x nTr x
%        nCh. The order dimention is chosen in a way that be more handy and
%        general
% 2      nBins: number of bins based on input binSize
%
% ------
% potential improvments:
% In case everything is based on sameple maybe good idea to have output of
% time axis 
% ------
% Code Info:
%   creation: 2013-05-04 by ShS -> shervin.safavi@gmail.com
%   modification: 2014-31-07 by ShS

%% data parameters 
[nCh, nTr] = size(spikeTrains);

% <<<<<<< HEAD
%% bin information
startT = times.startT;
stopT = times.stopT;
% binSize = 100;
binCenters =  startT + binSize / 2 :  binSize : stopT - binSize / 2;
nBin = numel(binCenters);
% =======
% %% bin structure
% % the initial and end point of interest; it can also be asked as an input
% startT = 0;
% stopT = 20000;
% % bins
% binCenters =  startT + binSize / 2 :  binSize : stopT - binSize / 2;
% nBins = numel(binCenters);

% %% spike trian
% if isempty(spikeTrain)
%     % considering the spike count for all bins zero
%     binSpikeCount = zeros(1, nBins);
% else
% % >>>>>>> origin/FRestimate

%% extracting firing rates
spikeCounts = zeros(nBin, nTr, nCh);
for iCh = 1 : nCh
    for iTr = 1 : nTr
        if isempty(spikeTrains{iCh, iTr})
            % considering the spike count for all bins zero
            spikeCounts(:, iTr, iCh) = zeros(nBin, 1);
        else            
            tmpRes = spikeTrains{iCh, iTr};
            nSpikes = numel(tmpRes);
            
            % makeing "nSpikes x 2" bin matrix (look 'bsxfun' in matlab document for more explnation)
            bins = bsxfun(@plus,[-binSize / 2 binSize / 2], binCenters');
            
            % ~ memory allocation
            spikeLogicalAsc = zeros(nSpikes, nBin);
            for binNum = 1 : nBin
                spikeLogicalAsc(:, binNum) = (bins(binNum, 1) <= tmpRes & tmpRes < bins(binNum, 2));
            end
            
            spikeCounts(:, iTr, iCh) = sum(spikeLogicalAsc, 1);
        end
    end
end