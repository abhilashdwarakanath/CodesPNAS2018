function [S, f] = lfpspectra(params,lfp,tag)

for i = 1:size(lfp,1)
    
    switch tag
        
        case 1 %on
            
            for tt=1:params.ntrials
                st1(:,tt) = (lfp{i,tt}(1:params.fulltime))* 2e-3;
            end
            
        case 0 %off
            
            for tt=1:params.ntrials
                st1(:,tt) = (lfp{i,tt}(params.fulltime+1:2*params.fulltime+30))* 2e-3;
            end
            
    end
    
   % Compute spectrum
   
   [S(:,i),f]=mtspectrumc(st1,params);
    
end

end