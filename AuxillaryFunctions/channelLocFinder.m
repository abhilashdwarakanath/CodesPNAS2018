function [MUArow MUAcol] = channelLocFinder(channelNumbers, channelsMap, MUAmap)
% [MUArow MUAcol] = channelLocFinder(channelNumbers, channelsMap, MUAmap)
% 
% ------
% Input:
% 
% Output:
% 
% ------
% Code Info:
%   creation: uknown 2013 or 2012 by Shervin Safavi (ShS) shervin.safavi@gmail.com
%   modification: 2014-08-03 by ShS



 nCh = numel(channelNumbers);

for i = 1 : nCh
    tmpLfpChIndex = channelsMap(channelNumbers(i) == channelsMap(:, 3), 2);
    [MUArow(i), MUAcol(i)] = find(MUAmap == tmpLfpChIndex);
%     mappedImage(MUAmap == tmpBoardIdex) = imgVec(channelNum);
end