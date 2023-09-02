function angle = angle_fromXaxis(p1, p2, varargin)
% angle = angle_fromXaxis(p1, p2, angleUnit)
% Simply computing the angle between line made by p1 and p2 and the
% x-axis. 
% Note that use this function only for 2D case, for general case, the
% general defenition of angle between two vector should be used.
%
% EXAMPLE:
% angle_fromXaxis([1 1], [2 2])
%
% ------
% Input:
% 1,2   p1,p2: 1*2 or 2*1 vectors
% (3)   angleUnit: string
%       define the unit and range the angle should be returened
%           'rad':  [-pi pi] 
%           'rad+': [0 2*pi] 
%           'deg':  [0 360]
% 
% Output:
% 1     angle: scalar
%
% ------
% see also atan2
% ------
% 
% potential improvments:
% (1) 
% ------
% Code Info:
%   creation: 2015-06-29 by ShS (shervin.safavi@gmail.com)
%   modification:
%       $ 201?-??-?? ?

%% Handle optional inputs (varargin): angleUnit
optionalVariables.angleUnit = []; 
defaultValues{1} = 'deg'; 
optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%%
p1_x = p1(1); p1_y = p1(2);
p2_x = p2(1); p2_y = p2(2);

angle_rad = atan2(p2_y - p1_y, p2_x - p1_x);


switch optionalVariables.angleUnit
    case 'deg' % [0 360]
        angle_deg = rad2deg(angle_rad);
        angle = angle_deg;
    case 'rad', angle = angle_rad; % [-pi pi] 
    case 'rad+' % [0 2*pi] 
        if angle_rad < 0, angle = angle_rad + 2*pi;
        else angle = angle_rad; end
    otherwise
end


