function [meanReplacedData binCenters] = mean_inBins(data, binSize, varargin)
% [meanTrend binCenters] = mean_inBins(data, binSize, locations)
% 
% ------
% Input:
% 1      i1
% 
% Output:
% 1      o1
%
% ------
% potential improvments:
% ------
% Code Info:
%   creation: 2014-08-17 by ShS -> shervin.safavi@gmail.com
%   modification: 

%% initialization

if nargin > 2
    locations = varargin{1};   
end

% assigning the diffault values
if ~exist('locations', 'var') 
    locations.start = min(data(:));
    locations.stop = max(data(:));
end

%% bin information
start = locations.start;
stop = locations.stop;
% binCenters =  start + binSize / 2 :  binSize : stop - binSize / 2;
binCenters =  start - binSize / 2 :  binSize : stop + 3*binSize / 2;
nBin = numel(binCenters);

% makeing "nSpikes x 2" bin matrix (look 'bsxfun' in matlab document for more explnation)
bins = bsxfun(@plus,[-binSize / 2 binSize / 2], binCenters');

%% find the means within bins and replace them with bin centers
meanReplacedData = data;
for iBin = 1 : nBin
    tmp_index = find(bins(iBin, 1) <= data & data < bins(iBin, 2));
    meanReplacedData(tmp_index) = binCenters(iBin);
end

