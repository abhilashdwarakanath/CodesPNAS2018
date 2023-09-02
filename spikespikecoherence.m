function [ssc,  sscohgram,sswavelet] = spikespikecoherence(params,cell_pair,spikeTimes,i,tag)

for r = 1:size(cell_pair{i},1)
    
    switch tag
        
        case 1 %on
            
            for tt=1:params.numtrials
                st1{1, tt} = spikeTimes{cell_pair{i}(r,1),tt}(spikeTimes{cell_pair{i}(r,1), tt}(:, 1)<=params.fulltime, 1);
                spikeTrain1(tt,:) = zeros(1,params.fulltime);
                if ~isempty(st1{1, tt})
                    spikeTrain1(tt,st1{1,tt}) = 1;
                end
                st2{1, tt} = spikeTimes{cell_pair{i}(r,2),tt}(spikeTimes{cell_pair{i}(r,2), tt}(:, 1)<=params.fulltime, 1);
                spikeTrain2(tt,:) = zeros(1,params.fulltime);
                if ~isempty(st2{1, tt})
                    spikeTrain2(tt,st2{1,tt}) = 1;
                end
                
                st1{1,tt} = st1{1,tt}*2e-3;
                st2{1,tt} = st2{1,tt}*2e-3;
            end
            
        case 0 %off
            
            for tt=1:params.numtrials
                st1{1, tt} = spikeTimes{cell_pair{i}(r,1),tt}(spikeTimes{cell_pair{i}(r,1), tt}(:, 1)>params.fulltime, 1)-params.fulltime;
                spikeTrain1(tt,:) = zeros(1,params.fulltime+30);
                if ~isempty(st1{1, tt})
                    spikeTrain1(tt,st1{1,tt}) = 1;
                end
                st2{1, tt} = spikeTimes{cell_pair{i}(r,2),tt}(spikeTimes{cell_pair{i}(r,2), tt}(:, 1)>params.fulltime, 1)-params.fulltime;
                spikeTrain2(tt,:) = zeros(1,params.fulltime+30);
                if ~isempty(st2{1, tt})
                    spikeTrain2(tt,st2{1,tt}) = 1;
                end
                
                st1{1,tt} = st1{1,tt}*2e-3;
                st2{1,tt} = st2{1,tt}*2e-3;
            end
            
    end
    
    % Compute coherence
    
    % Check for null activity cells
    
    emptytrls1 = cellfun(@isempty,st1);
    emptytrls2 = cellfun(@isempty,st2);
    
    null1 = sum(emptytrls1); null2 = sum(emptytrls2);
    
    if null1==200 && null2 == 200
        
        ssc{r} = NaN;
        sscohgram{r} = NaN;
        
    else
        
        % Convert to Chronux compatible form
        st1Chr = convert_spike_chronux(st1);
        st2Chr = convert_spike_chronux(st2);
        
        % Compute wavelet-coherogram
        for tr = 1:params.numtrials
            [sswcoh(:,:,tr),sswphase(:,:,tr),wf(1,:)] = wcoherence(spikeTrain1(tr,:),spikeTrain2(tr,:),500);
        end
        
        sswavelet{r}.coherogram = nanmean(sswcoh,3); sswavelet{r}.phase = nanmean(sswphase,3); sswavelet{r}.f = wf;
        
        % Compute coherence
        
        [ssc{r},~,~,~,~,freqs]=coherencypt(st1Chr.trial,st2Chr.trial,params);
        movingwin = [200e-3 20e-3];
        [C,phi,~,~,~,ctime]=cohgrampt(st1Chr.trial,st2Chr.trial,movingwin,params,1);
        sscohgram{r}.coherogram = C; sscohgram{r}.phase = phi; sscohgram{r}.ctime = ctime; sscohgram{r}.f = freqs;
        
        msg=sprintf('Valid pair : %d of %d\n',r,size(cell_pair{i},1));
        
        disp(msg)
    end
    
    %mssc = mean(ssc,2);
    
end

