function indices = find_zeroCrossing(v)
% indices = find_zeroCrossing(v)
% Find points in vector 'v' where sign between consecutive elemnts changes.
% i.e. points nearest to zero-crossings.
% In case you want to find the points where cross specific values you
% should substract that value from the vector v and then use the function.
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
% >> find_zeroCrossing(A)
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
% 1     v: input 1D vector which you need to find the zero-crossing points
%       
% Output: 
% 1     indices: indicies of zero-crossing points
%
% ------
% potential improvments:
% (1) find the value-crossing (not just zero), substract the value at begining and them
% do the same
% (2) zero crossing can be from + to - or - to + which potentially can be
% sperated and implemented in the inputs
% ------
% Code Info:
%   creation: 2015-01-16 but stolen from:
%   http://hips.seas.harvard.edu/content/count-zero-crossings-matlab
%   modification: 
%       $ 2015
%

% here you can find the complicated code :D
indices = find(diff(v > 0) ~= 0) + 1;
  
