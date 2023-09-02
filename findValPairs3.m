function cell_pairs = findValPairs3(distangs,corrmat,SUind)

for u = 1:size(distangs,2)
        
        [x y] = ind2sub(size(corrmat),distangs{u});
        
        clear cell_pair
        cell_pair = [];
        
        for tt=1:size(x,1);
            
            if ismember(x(tt,1),SUind)==1 && ismember(y(tt,1),SUind)==1
                %cell_pairs1 = sub2ind(size(corrmat),x(tt,1),y(tt,1));
                cell_pairs1 = [x(tt,1) y(tt,1)];
                cell_pair = [cell_pair;cell_pairs1];
            end
        end
        
        
        cell_pairs{u} = cell_pair;
end