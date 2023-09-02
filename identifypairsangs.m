function distangs = identifypairsangs(SU_angles,SU_PWdistances)

%% Identify pairs with same separations and save the coordinates of all pairs in a particular angsance bin

ang_bin_bounds = 0:60:360;
%ang_bin_bounds = 0:pi/3:2*pi;
dist_bin_bounds = 0.25:0.5:4.5;

angs = cell(1,length(ang_bin_bounds)-1);
for k = 1:length(ang_bin_bounds)-1
    
    angs{k}=find(SU_angles >= deg2rad(ang_bin_bounds(k)) & SU_angles < deg2rad(ang_bin_bounds(k+1)));
    
end

dist = cell(1,length(dist_bin_bounds)-1);
for k = 1:length(dist_bin_bounds)-1
    
    dist{k}=find(SU_PWdistances >= (dist_bin_bounds(k)) & SU_PWdistances < (dist_bin_bounds(k+1)));
    
end

distangs = cell(length(dist),length(angs));

for i = 1:length(dist)
    for j = 1:length(angs)
        distangs{i,j} = intersect(dist{i},angs{j});
    end
end