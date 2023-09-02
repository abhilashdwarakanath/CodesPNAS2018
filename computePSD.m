function [meanpsd,f] = computePSDlfp(data,decimFac,newLength,numPieces,trlInSamps,dps_seq,lambda)

decsig = decimate(im2double(data,decimFac));
    decsig = decsig(1:newLength);
    decsig = reshape(decsig,trlInSamps,numPieces);
    
    [p,f] = pmtm(decsig,dps_seq,lambda',[],500,'eigen');
    meanpsd = mean(10*log10(p),2);
end