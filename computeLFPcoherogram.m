function [wcoh,fcoher] = computeLFPcoherogram(unit1,unit2)

%% Reshape both signals

%decsig1 = unit1(1:newLength);
%decsig2 = unit2(1:newLength);
%decsig1 = reshape(decsig1,trlInSamps,numPieces);
%decsig2 = reshape(decsig2,trlInSamps,numPieces);

decsig1 = unit1;
decsig2 = unit2;

%% Compute coherogram

% I chose these "scales" to get upto 250Hz. How do I justify this?

a0 = 2^(1/12);
scales = a0.^linspace(exp(log(0.000005*2)),exp(log(64)),250);

for pieces = 1:size(decsig1,2)
    [wcoh(:,:,pieces)] = wcoher(decsig1(:,pieces),decsig2(:,pieces),scales,'cmor1-0.1');
    %[wcoh(:,:,pieces),~,fcoher] = wcoherence(decsig1(:,pieces),decsig2(:,pieces),500);
end

%% Compute means

wcoh = nanmean(wcoh,3);

end
