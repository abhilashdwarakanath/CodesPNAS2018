function PWdistances = generateUtahPWdistances(integratedTabel)
% PWdistances = generateUtahPWdistances(integratedTabel);
% ------
% Input:
% 1      integratedTabel: This table should contain the channel/MU/SU 
%number and the its location in the array. integratedTabel is nCh/nSU/nMU 
%x 2 matrix with this columnar organization:
%        locationRowInArray #     locationRowInArray #
% Output:
% 1      PWdistances: All pair-wise distances
%
% ------
% potential improvments:
% (1) ask the channel distance unit as an  input
% (2) check the size of input it should be only 2 columns
% ------
% Code Info:
%   creation: 2014-08-17 by ShS -> shervin.safavi@gmail.com
%   modification:
%
%% initialization
chDistanceUnit = .4; % Center To Center Distance of two neighbour 
%channels in mm
nCh = size(integratedTabel, 1);
% [nArrayRows nArrayColumns]= size(arrayMap);

PWdistances = zeros(nCh, nCh);

%% compute euclidean distance
% you don't need to compute for all pairs. As
% the distance is symetric you can only compute it for up or down triabgle
% of the sqaure matrix; but here is for the general case
for iCh = 1 : nCh
    for jCh = iCh+1 : nCh
        % compute the pair-wise distances
        PWdistances(iCh, jCh) = ...
            norm([integratedTabel(iCh, 1) integratedTabel(iCh, 2)]*chDistanceUnit ...
            - [integratedTabel(jCh, 1) integratedTabel(jCh, 2)]*chDistanceUnit);
        PWdistances(jCh, iCh) = PWdistances(iCh, jCh);
    end
end
