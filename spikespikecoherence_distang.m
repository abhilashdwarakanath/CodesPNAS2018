function [ssc, freqs] = spikespikecoherence_distang(params,cell_pair,spikeTimes,tag)

for r = 1:size(cell_pair,1)
    
    switch tag
        
        case 1 %on
            
            for tt=1:params.numtrials
                st1{1, tt} = spikeTimes{cell_pair(r,1),tt}(spikeTimes{cell_pair(r,1), tt}(:, 1)<=params.fulltime, 1)* 2e-3;
                st2{1, tt} = spikeTimes{cell_pair(r,2),tt}(spikeTimes{cell_pair(r,2), tt}(:, 1)<=params.fulltime, 1)* 2e-3;
            end
            
        case 0 %off
            
            for tt=1:params.numtrials
                st1{1, tt} = spikeTimes{cell_pair(r,1),tt}(spikeTimes{cell_pair(r,1), tt}(:, 1)>params.fulltime, 1) * 2e-3;
                st2{1, tt} = spikeTimes{cell_pair(r,2),tt}(spikeTimes{cell_pair(r,2), tt}(:, 1)>params.fulltime, 1) * 2e-3;
            end
            
    end
    
    % Convert to Chronux compatible form
    st1Chr = convert_spike_chronux(st1);
    st2Chr = convert_spike_chronux(st2);
    
    % Compute coherence
    %addAttachedFiles('local', ssc)
    [ssc(:,r),~,~,~,~,freqs]=coherencypt(st1Chr.trial,st2Chr.trial,params);
     msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair,1));
    disp(msg)
end
end

