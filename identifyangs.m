function angs = identifyangs(SU_angles)

%% Identify pairs with same separations and save the coordinates of all pairs in a particular angsance bin

ang_bin_bounds = 0:60:360;

angs = cell(1,length(ang_bin_bounds)-1);
for k = 1:length(ang_bin_bounds)-1
    
    angs{k}=find(SU_angles >= deg2rad(ang_bin_bounds(k)) & SU_angles < deg2rad(ang_bin_bounds(k+1)));
    
end
end