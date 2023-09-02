function indices = find_valueCrossing(timeSeries, varargin)
% indices = find_valueCrossing(timeSeries, value)
% Find points in 'timeSeries' where sign [relative to 'value'] between consecutive elemnts changes.
% i.e. points nearest to value-crossings.
% 
% example:
% >> A = randn([10 1])
%
% A =
%
%  -0.49840598306643
%   1.04975509964655
%  -1.67055867973620
%  -2.01437026154355
%   0.98661592496732
%  -0.06048256273708
%   1.19294080740269
%   2.68558025885591
%   0.85373360483580
%   1.00554850567375
%
% >> find_valueCrossing(A, 0)
%
% ans =
%
%     2
%     3
%     5
%     6
%     7
% 
% ------
% Input:
% 1     timeSeries: 1D vector
%       input 1D vector which you need to find the zero-crossing points
% (2)   Value: scalar 
%       the crossing point (scalar)
%       default: 0
%       
% Output: 
% 1     indices: indicies of value-crossing points
%
% ------
% potential improvments:
% (1) zero crossing can be from + to - or - to + which potentially can be
% sperated and implemented in the inputs
% ------
% Code Info:
%   creation: 2015-01-25 by ShS -> shervin.safavi@gmail.com 
%   somehow is stolen from:
%   http://hips.seas.harvard.edu/content/count-zero-crossings-matlab
%   modification: 
%       $ 2015
%

%% handle optional inputs (varargin): value
% optional variables:
optionalVariables.value = [];
% default values:
defaultValues{1} = 0;

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%% here you can find the complicated code for finding the crossing points :D
timeSeries_valueSubstracted = timeSeries - optionalVariables.value;
indices = find(diff(timeSeries_valueSubstracted  > 0) ~= 0) + 1;
  
