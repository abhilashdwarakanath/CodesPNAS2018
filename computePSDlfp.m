function [meanpsd,f] = computePSDlfp(data,newLength,numPieces,trlInSamps,dps_seq,lambda,targetFreq,flag)

decsig = data(1:newLength);
decsig = reshape(decsig,trlInSamps,numPieces);

switch flag
    
    case 1
        
        [p,f] = pmtm(decsig,dps_seq,lambda,[],targetFreq);
        %meanpsd = mean(10*log10(p),2);
        p = sum(p,2);
        meanpsd=tsmovavg(p','e',10);
        
    case 2
        
        win = barthannwin(size(decsig,1));
        [p,f] = periodogram(decsig,win,'onesided',500,targetFreq);
        %meanpsd = mean(10*log10(p),2);
        p = sum(p,2);
        meanpsd=tsmovavg(p','e',10);
        
end
end