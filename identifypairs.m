function dist = identifypairs(SU_PWdistances)

%% Identify pairs with same separations and save the coordinates of all pairs in a particular distance bin

dist_bin_bounds = 0.25:0.5:4.5;

dist = cell(1,length(dist_bin_bounds)-1);
for k = 1:length(dist_bin_bounds)-1
    
    dist{k}=find(SU_PWdistances >= (dist_bin_bounds(k)) & SU_PWdistances < (dist_bin_bounds(k+1)));
    
end
