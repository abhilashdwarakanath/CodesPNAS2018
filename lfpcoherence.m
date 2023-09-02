function [cohon, cohoff] = lfpcoherence(params,lfp,i)

cohon = zeros(params.units,params.units,length(params.oscillations)-1);
cohoff = zeros(params.units,params.units,length(params.oscillations)-1);

for j = 1:params.units %nCells
    
    for k = 1:params.units %nCells
        fprintf('Computing coherence for Dataset : %d, pairs %d and %d\n',i,j,k)
        tic;
        xon = zeros(params.trialtime,params.ntrials);
        yon = zeros(params.trialtime,params.ntrials);
        xoff = zeros(params.trialtime+30,params.ntrials);
        yoff = zeros(params.trialtime+30,params.ntrials);
        for l = 1:params.ntrials %nTrials
            xon(:,l) = lfp{j,l}(1:params.trialtime); % Stim ON 10s First Unit
            yon(:,l) = lfp{k,l}(1:params.trialtime); % Stim ON 10s Second Unit
            xoff(:,l) = lfp{j,l}(params.trialtime+1:(params.trialtime*2)+30); % Stim OFF 10s First Unit
            yoff(:,l) = lfp{k,l}(params.trialtime+1:(params.trialtime*2)+30); % Stim ON 10s Second Unit
        end
        % Compute pair-wise coherence measure in 5 bands - Delta, Theta,
        % Beta and Gamma. Compute the "energy" of coherence by taking the
        % mean
        for bands = 1:length(params.oscillations)-1
            params.fpass = [params.oscillations(bands) params.oscillations(bands+1)];
            con=coherencyc(xon,yon,params);
            coff=coherencyc(xoff,yoff,params);
            con(con==0) = nan;
            coff(coff==0) = nan;
            cenergyon = nanmean(con);
            cenergyoff = nanmean(coff);
            cohon(j,k,bands) = cenergyon;
            cohoff(j,k,bands) = cenergyoff;
        end
        toc;
    end
    
end