function [C,f] = computeLFPcoherence(unit1,unit2,targetFreq,tapers,newLength,trlInSamps,numPieces,method)

%% Reshape both signals

newLength = 15000;
numPieces = 3;

decsig1 = unit1(1:newLength);
decsig2 = unit2(1:newLength);
decsig1 = reshape(decsig1,trlInSamps,numPieces);
decsig2 = reshape(decsig2,trlInSamps,numPieces);

%decsig1 = unit1;
%decsig2 = unit2;

%% Compute coherence. Currently I am using mscohere, but later on I will implement the Chronux function

switch method
    
    case 1 % Hanning window
        
        [C,f] = mscohere(decsig1,decsig2,barthannwin(size(decsig1,1)),[],[],targetFreq);
        
        C = sum(C,2);
        C = tsmovavg(C','e',10); % 10 points amounts to 1 Hz
        
    case 2 % Multitaper
        
        for pieces = 1:size(decsig1,2)
            
            [C(:,pieces), f] = cmtm(decsig1(:,pieces),decsig2(:,pieces),tapers,targetFreq);
            
        end
        
        C = abs(C);
        C = nanmean(C,2);
        
end

end