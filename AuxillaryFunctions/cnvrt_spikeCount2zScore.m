function zScores = cnvrt_spikeCount2zScore(spikeCounts)
% zScores = cnvrt_spikeCount2zScore(spikeCounts)
% It's just a simple helper function, that covert spike count to z-score
% for each condition
%       z = \frac{spikeCount - mean}{std}
% 
% All these for one condition
% 
% The same function can be applied on firing rates
% 
% ------
% Input:
% 1      spikeCounts: matrix of 'spikeCounts' which should have this
%        dimention:
%           nBins * nTr * nUnit
%           (nUnit can be number of channels, single unit etc)
% 
% Output:
% 1      z-score values with the same dimention as input
%
% ------
% potential improvments:
% (1) check the input dimention
% ------
% Code Info:
%   creation: 2014-10-29 by ShS -> shervin.safavi@gmail.com
%   modification:
%       $ 

%% extract the dimentions
[n.Bins, n.Tr, n.Unit] = size(spikeCounts);
zScores = nan(size(spikeCounts));

%% converting to z-score
for iBin = 1 : n.Bins
    for iUnit = 1 : n.Unit
        zScores(iBin, :, iUnit) = cnvrt_x2zScore(spikeCounts(iBin, :, iUnit));
    end 
end