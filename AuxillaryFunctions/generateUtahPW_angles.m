function PWanglesFromXaxis = generateUtahPW_angles(integratedTabel)
%% UTAH ARRAY FUNCTION 
% PWanglesFromXaxis = generateUtahPW_angles(integratedTabel)
% Will compute the angle between line made by two channels/units and the
% x-axis. This x-axis is already defined by the coordinate (row-column
% index) you address the channels/units.
% 
% You are assmung that the center the coordinate system is one
% the top-right corner. 
% y+: donward
% x+: rightward
% The coordinate system has been transformed in the function:
% new y+: upward
% new origin: buttom-left
% 
% EXAMPLE:
% generateUtahPW_angles(utahMaps(iDataset).MU.integratedTabel(:, [3 4]));
%
% ------
% Input:
% 1     integratedTabel: nCh/Unit * 2 matrix
%       This table should contain the  channel/MU/SU number and their location in the array  
%       The one which is assumed is this code is supposed to be generated by 'generatUtahMaps' function which is nCh x (~)5 matrix with this columnar organization:
%       board#    LFP#    spk#    locationRowInArray #     locationRowInArray #
%       The function ONLY needs the 2 columns of the integratedTabel that
%       indicate the location of the channel/unit in the array
% 
% Output:
% 1     PWanglesFromXaxis: 
%       All pair-wise distances
%
% ------
% see also generateUtahPWdistances generatUtahMaps2
% ------
% potential improvments:
% 
% ------
% Code Info:
%   creation: 2015-06-30 by ShS (shervin.safavi@gmail.com)
%   modification:
%       $ 201?-0?-?? ?

%% Some Initialization
% chDistanceUnit = .4; % Center To Center Distance of two neighbour channels in mm
nCh = size(integratedTabel, 1);
PWanglesFromXaxis = nan(nCh, nCh);

nCh_Y = 10;
% *** might need to be changed, it's a hardcoded value

%% Angle from X-axis
for iCh = 1 : nCh
    for jCh = 1 : nCh
        % compute the angles        
        Xi = integratedTabel(iCh, 2);         
        Yi = -1 * integratedTabel(iCh, 1) + nCh_Y+1; % Y should be transformed to new coordiante system
        
        Xj = integratedTabel(jCh, 2); 
        Yj = -1 * integratedTabel(jCh, 1) + nCh_Y+1; % Y should be transformed to new coordiante system
                
%         PWanglesFromXaxis(iCh, jCh) = atan2(Yj-Yi, Xj-Xi); 
        PWanglesFromXaxis(iCh, jCh) = angle_fromXaxis([Xi Yi], [Xj Yj], 'rad+');    
    end
end