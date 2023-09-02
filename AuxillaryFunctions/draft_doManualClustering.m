function cluIndex=draft_doManualClustering(HDdata, Fvec)
% draft_doManualClustering(HDdata, Fvec)

 cluIndex = recordClickTrace(Fvec(1,:), Fvec(2,:));
 figure(102)
 clf
%  HDdata_nrmlz = bsxfun(@rdivide, HDdata,  max(HDdata));
%  plot(HDdata_nrmlz(:, cluIndex))
%plot(HDdata(:, cluIndex),'k')
%hold all
%plot(nanmean(HDdata(:, cluIndex), 2), 'LineWidth',5,'Color',[1 0 1]);