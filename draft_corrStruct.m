% a draft for extractong corrolation structres
% exploring in correlation structures
%% plotting firing rates
FR = mean(allSpikeCounts, 2);
FR = squeeze(FR);
subplot_wUtahMap(FR, 'linearFR', 'Firing Rates');
plot(FR(:, 94), 'k');

%% play w/ matlab corr function

%% see some corrolation in the array
% accumulating the trials 
allSpikeCounts_allTrAcum = reshape(allSpikeCounts, size(allSpikeCounts,1)*nTr, nCh);

% a chosen column of array
chVertical = MUAmap(:, 7);
chHorizental = MUAmap(3, :);
chDiag = diag(MUAmap); chDiag = chDiag(~isnan(chDiag)); 

% compute the corrolation in both cases
% vertical
for iCh = 2 : numel(chVertical)
   R = corrcoef(allSpikeCounts_allTrAcum(:, chVertical(1)), allSpikeCounts_allTrAcum(:, chVertical(iCh))); 
   corrVer(iCh-1) = R(2,1);
end

% vertical
for iCh = 2 : numel(chHorizental)
   R = corrcoef(allSpikeCounts(:,:, chHorizental(1)), allSpikeCounts(:,:, chHorizental(iCh))); 
   corrHor(iCh-1) = R(2,1);
end

% diag
for iCh = 2 : numel(chDiag)
   R = corrcoef(allSpikeCounts(:,:, chHorizental(1)), allSpikeCounts(:,:, chHorizental(iCh))); 
   corrDiag(iCh-1) = R(2,1);
end

subplot 311
plot(corrVer, 'k')
subplot 312
plot(corrHor, 'k')
subplot 313
plot(corrDiag, 'k')

%%  some example of 2D corr
% let's pick one channel in the middle of the array as the reference and
% see how's corrolated with others

% all channels
[MUArow MUAcol] = channelLocFinder((1:nCh), channelsMap, MUAmap);

refCh = 1;

corrStruct = zeros(10);
for iCh = 1 : nCh
    %     if ~isnan(MUAmap(MUAmap == iCh))
    %% corr structure
    R = corrcoef(allSpikeCounts(:,:, refCh), allSpikeCounts(:,:, iCh));
    corrStruct(MUArow(iCh), MUAcol(iCh)) = R(1,2);
    %% gradient of corr structure
    [csgX,csgY] = gradient(corrStruct);
    abs_csg = csgX.^2 + csgY.^2;
    %% distance stuff
    chDistanceFromRef(iCh, 1) = norm([MUArow(iCh) MUAcol(iCh)] - [MUArow(refCh) MUAcol(refCh)]);
    allCorrVal(iCh, 1) = R(1,2);
    %     end
end
% find the unique distances
uniqueDistances = unique(chDistanceFromRef);

for iDist = 1 : numel(uniqueDistances)
    [tmp_r ~] = find(chDistanceFromRef == uniqueDistances(iDist));
    clear tmp__corr_Func_distance
    tmp__corr_Func_distance(:, 1) = chDistanceFromRef(tmp_r); 
    tmp__corr_Func_distance(:, 2) = allCorrVal(tmp_r);
    corr_Func_distance(iDist, :) = mean(tmp__corr_Func_distance, 1);
end

% [chDistanceSorted, ndx_ascendingDistance] = sort(chDistanceFromRef);
% corr_Func_distance(:, 1) = chDistanceSorted;
% corr_Func_distance(:, 2) = allCorrVal(ndx_ascendingDistance);
corr_Func_distance(1, 2) = NaN;

% [r c] = find(MUAmap == refCh);
corrStruct(MUArow(refCh), MUAcol(refCh)) = 0;
clf
subplot 221
hold all
imagesc(corrStruct); axis equal; axis off; colormap gray
plot(MUAcol(refCh), MUArow(refCh), '+', 'color', 'r')
subplot 222
hold all
imagesc(abs_csg); axis equal; axis off; colormap gray
quiver(1:10, 1:10, csgX, csgY, 'LineWidth', 2, 'color', 'g')
plot(MUAcol(refCh), MUArow(refCh), 'o', 'color', 'r')
subplot(2,2, [3 4])
hold all
plot(corr_Func_distance(:, 2), 'color', 'b')
% subplot(3,2, [5 6])
% plot(sign(diff(corr_Func_distance(:, 2))), '.')
plot(smooth(corr_Func_distance(:, 2)), 'color', 'k', 'LineWidth', 2)
hold off

%% [VOID] integrative analysis for all the channels 
% all channels
[MUArow MUAcol] = channelLocFinder((1:nCh), channelsMap, MUAmap);

for refCh = 1 : nCh    
%     corrStruct = zeros(10);
    for iCh = 1 : nCh
        %% corr structure
        R = corrcoef(allSpikeCounts(:,:, refCh), allSpikeCounts(:,:, iCh));
        corrStruct(MUArow(iCh), MUAcol(iCh)) = R(1,2);
%         %% gradient of corr structure
%         [csgX,csgY] = gradient(corrStruct);
%         abs_csg = csgX.^2 + csgY.^2;
        %% distance stuff
        chDistanceFromRef(iCh, 1) = norm([MUArow(iCh) MUAcol(iCh)] - [MUArow(refCh) MUAcol(refCh)]);
        allCorrVal(iCh, 1) = R(1,2);
    end
    % find the unique distances
    uniqueDistances = unique(chDistanceFromRef);
    
    for iDist = 1 : numel(uniqueDistances)
        [tmp_r ~] = find(chDistanceFromRef == uniqueDistances(iDist));
        clear tmp__corr_Func_distance
        tmp__corr_Func_distance(:, 1) = chDistanceFromRef(tmp_r);
        tmp__corr_Func_distance(:, 2) = allCorrVal(tmp_r);
        corr_Func_distance(iDist, :) = mean(tmp__corr_Func_distance, 1);
    end
    
    % [chDistanceSorted, ndx_ascendingDistance] = sort(chDistanceFromRef);
    % corr_Func_distance(:, 1) = chDistanceSorted;
    % corr_Func_distance(:, 2) = allCorrVal(ndx_ascendingDistance);
    corr_Func_distance(1, 2) = NaN;
    
end

%% correlation structure of all channels together
allSpikeCounts_allTrAcum = reshape(allSpikeCounts, size(allSpikeCounts,1)*nTr, nCh);
R = corr(allSpikeCounts_allTrAcum);
allChPWdistances = generateUtahPWdistances(utahMaps.MUA.integratedTabel(:, [3 4 5]));

% sort the distance and their corresponding corrlations
% find the unique distances
uniqueDistances = unique(allChPWdistances);

for iDist = 1 : numel(uniqueDistances)
    tmp_index = find(allChPWdistances == uniqueDistances(iDist));
    clear tmp__corr_Func_distance
    tmp__corr_Func_distance(:, 1) = allChPWdistances(tmp_index);
    tmp__corr_Func_distance(:, 2) = R(tmp_index);
    corr_Func_distance(iDist, :) = mean(tmp__corr_Func_distance, 1);
end

hold all
plot(corr_Func_distance(:, 2), 'color', 'b')
plot(smooth(corr_Func_distance(:, 2)), 'color', 'k', 'LineWidth', 2)

%% divide to bins 
allSpikeCounts_allTrAcum = reshape(allSpikeCounts, size(allSpikeCounts,1)*nTr, nCh);
R = corr(allSpikeCounts_allTrAcum);
allChPWdistances = generateUtahPWdistances(utahMaps.MUA.integratedTabel(:, [3 4 5]));

% put the distances in bins with 0.4 mm steps
[meanReplacedDistances distance_binCenters] = mean_inBins(allChPWdistances, 0.05);

uniqueDistances = unique(meanReplacedDistances);
clear corr_Func_distance
for iDist = 1 : numel(uniqueDistances)
    tmp_index = find(meanReplacedDistances == uniqueDistances(iDist));
    clear tmp__corr_Func_distance
    tmp__corr_Func_distance(:, 1) = meanReplacedDistances(tmp_index);
    tmp__corr_Func_distance(:, 2) = R(tmp_index);
    corr_Func_distance(iDist, :) = mean(tmp__corr_Func_distance, 1);
end

figure
hold all
% plot(corr_Func_distance(:, 1), corr_Func_distance(:, 2), 'color', 'b')
plot(corr_Func_distance(:, 2), 'color', 'b')
plot(smooth(corr_Func_distance(:, 2)), 'color', 'k', 'LineWidth', 2)



%% look at all pairwise corrolation
R = corr(allSpikeCounts_allTrAcum);
R(R == 1) = 0;
imagesc(R); axis square; axis off
