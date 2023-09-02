function [unitRow unitCol] = channelLocFinder2(unitNumbers, unitMap)
% [MUArow MUAcol] = channelLocFinder(channelNumbers, channelsMap, MUAmap)
% 
% EXAMPLE:
%
% 
% ------
% Input:
% 
% Output:
% 
% ------
% Code Info:
%   creation: 2014-11-11 by ShS -> shervin.safavi@gmail.com
%   modification: 
%       $ 201? 



 nUnit = numel(unitNumbers);

for i = 1 : nUnit
    [unitRow(i), unitCol(i)] = find(unitMap == unitNumbers(i));
end