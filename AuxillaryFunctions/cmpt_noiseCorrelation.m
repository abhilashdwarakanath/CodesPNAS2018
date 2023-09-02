function correlationMatrix = cmpt_noiseCorrelation(data)
% correlationMatrix = cmpt_noiseCorrelation(data, normalizationFlag)
% 
% 
% example:
% 
% 
% ------
% Input:
% 1     zScores: nObservation * nVariable
%       a 2D array in which each row contain the observation/measurment of all variables and
%       each column contain all observations/measurments belong to
%       one variable.
% (2)   normalization: if all the observations/measurments are
%       normalized (tranfered to z-scores), the flag would be 1; then
%       mean won't be substracted and divided by vaiance
%       
% Output: 
% 1     correlationMatrix: nVariable * nVariable
%       the linear correlation for all pairs                   
%
% ------
% potential improvments:
% (1) what MATLAB 'corr' function does have extra
% (2) add an option if the data is not normalized 
% ------
% Code Info:
%   creation: 2014-10-30 by ShS -> shervin.safavi@gmail.com
%   modification: 
%       $ 2014-11-1 divison by number of observations was misssing

nObservation = size(data, 1);

zScores = data;
correlationMatrix = (1/nObservation)* zScores' * zScores;