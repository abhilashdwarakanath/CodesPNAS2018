%clear all
clc
%close all

% This program analyses various things from the GPFA stuff. Add proper
% description later.

%% load results and data

%load GPFA_2LV_100msbins_5Hz_resid_Results.mat
%load('spikeData20150813_UATE2010');

ndatasets = length(dataset);
dBins = 0.5:0.5:4;
aBins = deg2rad(30:60:330);
binCentersON.r = dBins;
binCentersON.theta = aBins-deg2rad(10);

binCentersOFF.r = dBins;
binCentersOFF.theta = aBins+deg2rad(10);

%% Compute distance bins etc

for currset = 1:params.ndatasets
    
    nDataSet = currset;
    msg=sprintf('Data set %d of %d starts...\n',currset,params.ndatasets);
    disp(msg);
    
    tic;
    
    % Extract spike times, unit distances
    
    spikeTime = jMU(1,currset).spikesBySampleIndex;
    nUnits = length(spikeTime);
    
    MU_PWdistances = tril(jMU(currset).PWdistances); % Remove doubles
    MU_PWangles = tril(jMU(currset).PWanglesFromXaxis); % Remove doubles
    
    % Identify pairs in distance bins
    
    dist1 = identifypairsangs(MU_PWangles,MU_PWdistances);
    
    for kk = 1:length(spikeTime)
        for j = 1:params.ntrials
            spikeTimes{kk,j}=spikeTime{kk}{j}; % spike times
        end
    end
    MUcontam = utahMaps(1,currset).MU.integratedTabel(:,5);
    SUind = find(MUcontam<=1000000);
    dist{currset} = findValPairs2(params,dist1,spikeTimes,SUind);
    
     % Compute PSTH
    psthON = zeros(nUnits,length(edgesON)-1,params.ntrials);
    psthOFF = zeros(nUnits,length(edgesOFF)-1,params.ntrials);
    respsthON = zeros(nUnits,length(edgesON)-1,params.ntrials);
    respsthOFF = zeros(nUnits,length(edgesOFF)-1,params.ntrials);
    
    for i = 1:nUnits
        for j = 1:params.ntrials
            spkson = spikeTime{i}{j}(spikeTime{i}{j}<=params.on);
            spksoff = spikeTime{i}{j}(spikeTime{i}{j}>params.on);
            if ~isempty(spkson)
                r1 = (histc(spkson,edgesON)/persec);
            else r1 = 1e-4*randn(1,length(edgesON))+0;
            end
            r1 = r1(1:end-1);
            if ~isempty(spksoff)
                r2 = (histc(spksoff,edgesOFF)/persec);
            else r2 = 1e-4*randn(1,length(edgesOFF))+0;
            end
            r2 = r2(1:end-1);
            
            mu1 = mean(r1);
            mu2 = mean(r2);
            psthON(i,:,j) = r1;
            psthOFF(i,:,j) = r2;
            respsthON(i,:,j) = r1-mu1;
            respsthOFF(i,:,j) = r2-mu2;
        end
    end
    
    % Compute naive factor loadings and noise
    
    rpon = reshape(respsthON,nUnits,(length(edgesON)-1)*params.ntrials);
    rpoff = reshape(respsthOFF,nUnits,(length(edgesOFF)-1)*params.ntrials);
    
    Qon = cov(rpon');
    Qoff = cov(rpoff');
    [naiveCon{currset}, Lambdaon] = eigs(Qon, 2);
    [naiveCoff{currset}, Lambdaoff] = eigs(Qoff, 2);
    
    naiveRon{currset} = diag(diag(Qon - naiveCon{currset} * Lambdaon * naiveCon{currset}'));
    naiveRoff{currset} = diag(diag(Qoff - naiveCoff{currset} * Lambdaoff * naiveCoff{currset}'));
    
    for i = 1:params.ndatasets
        map{i} = utahMaps(i).MU.map;
    end
end

%% Correlation structure

corrstrucdistON = zeros(params.lv,length(dBins),length(aBins));
corrstrucdistOFF = zeros(params.lv,length(dBins),length(aBins));
rem_corrstrucdistON = zeros(1,length(dBins),length(aBins));
rem_corrstrucdistOFF = zeros(1,length(dBins),length(aBins));

for currset = 1:ndatasets
    for nlv = 1:params.lv
        C_on{nlv} = dataset(currset).modelON.C(:,nlv)*dataset(currset).modelON.C(:,nlv)';
        C_off{nlv} = dataset(currset).modelOFF.C(:,nlv)*dataset(currset).modelOFF.C(:,nlv)';
    end
    
    residuals_on = dataset(currset).modelON.R;
    residuals_off = dataset(currset).modelOFF.R;
    rem_on = dataset(currset).modelON.rem;
    rem_off = dataset(currset).modelOFF.rem;
    
    % Compute correlation coefficients along distance
    
    corr_dist_on = zeros(nlv,size(dBins,2),size(aBins,2));
    corr_dist_off = zeros(nlv,size(dBins,2),size(aBins,2));
    rem_corr_dist_on = zeros(1,size(dBins,2),size(aBins,2));
    rem_corr_dist_off = zeros(1,size(dBins,2),size(aBins,2));
    
    for nBins = 1:size(dBins,2)
        for nAngs = 1:size(aBins,2)
            for nPairs = 1:size(dist{currset}{nBins,nAngs})
                for nLV = 1:params.lv
                    if ~isempty(dist{currset}{nBins,nAngs})
                        
                        [x,y] = ind2sub(size(C_on{nLV}),dist{currset}{nBins,nAngs}(nPairs));
                        
                        corr_dist_on(nLV,nBins,nAngs) = corr_dist_on(nLV,nBins,nAngs)+(C_on{nLV}(x,y)/sqrt(((C_on{nLV}(x,x))^2+residuals_on(x,x))*((C_on{nLV}(y,y))^2+residuals_on(y,y))));
                        corr_dist_off(nLV,nBins,nAngs) = corr_dist_off(nLV,nBins,nAngs)+(C_off{nLV}(x,y)/sqrt(((C_off{nLV}(x,x))^2+residuals_off(x,x))*((C_off{nLV}(y,y))^2+residuals_off(y,y))));
                        
                        %corr_dist_on(nLV,nBins,nAngs) = corr_dist_on(nLV,nBins,nAngs)+(C_on{nLV}(x,y)/sqrt(((C_on{nLV}(x,x))+residuals_on(x,x))^2*((C_on{nLV}(y,y))+residuals_on(y,y))^2));
                        %corr_dist_off(nLV,nBins,nAngs) = corr_dist_off(nLV,nBins,nAngs)+(C_off{nLV}(x,y)/sqrt(((C_off{nLV}(x,x))+residuals_off(x,x))^2*((C_off{nLV}(y,y))+residuals_off(y,y))^2));
                    end
                    
                end
                
                rem_corr_dist_on(1,nBins,nAngs) = rem_corr_dist_on(1,nBins,nAngs)+(rem_on(x,y)/sqrt(((rem_on(x,x))^2)*((rem_on(y,y))^2)));
                rem_corr_dist_off(1,nBins,nAngs) = rem_corr_dist_off(1,nBins,nAngs)+(rem_off(x,y)/sqrt(((rem_off(x,x))^2)*((rem_off(y,y))^2)));
                corr_dist_on(nLV,nBins,nAngs) = corr_dist_on(nLV,nBins)/length(dist{currset}{nBins,nAngs});
                corr_dist_off(nLV,nBins,nAngs) = corr_dist_off(nLV,nBins)/length(dist{currset}{nBins,nAngs});
                
            end
            
        end
        
        rem_corr_dist_on(1,nBins,nAngs) = rem_corr_dist_on(1,nBins,nAngs)/length(dist{currset}{nBins,nAngs});
        rem_corr_dist_off(1,nBins,nAngs) = rem_corr_dist_off(1,nBins,nAngs)/length(dist{currset}{nBins,nAngs});
        
    end
    
    dataset(currset).ON.corrstruc = corr_dist_on;
    dataset(currset).OFF.corrstruc = corr_dist_off;
    dataset(currset).ON.rem_corrstruc = rem_corr_dist_on;
    dataset(currset).OFF.rem_corrstruc = rem_corr_dist_off;
    
    corrstrucdistON = corrstrucdistON+corr_dist_on;
    corrstrucdistOFF = corrstrucdistOFF+corr_dist_off;
    rem_corrstrucdistON = rem_corrstrucdistON+rem_corr_dist_on;
    rem_corrstrucdistOFF = rem_corrstrucdistOFF+rem_corr_dist_off;
    
end

corrstrucdistON = corrstrucdistON/ndatasets;
corrstrucdistOFF = corrstrucdistOFF/ndatasets;
rem_corrstrucdistON = rem_corrstrucdistON/ndatasets;
rem_corrstrucdistOFF = rem_corrstrucdistOFF/ndatasets;

%% Spatial structures

% Plot correlation along distance, ON vs OFF, all latent variables

lw = [2 1];

% Compute SEM vector along distance
for nlv = 1:size(corrstrucdistON,1)
    errON(nlv,:) = sem(mean(corrstrucdistON(nlv,:,:),3))*ones(1,length(dBins));
    errOFF(nlv,:) = sem(mean(corrstrucdistOFF(nlv,:,:),3))*ones(1,length(dBins));
end

remerrON = sem(mean(rem_corrstrucdistON(1,:,:),3))*ones(1,length(dBins));
remerrOFF = sem(mean(rem_corrstrucdistOFF(1,:,:),3))*ones(1,length(dBins));

% Plot
figure(1)
for nlv = 1:params.lv
    subplot(2,2,nlv)
    errorbar(dBins,mean(corrstrucdistON(nlv,:,:),3),errON(nlv,:),'b','LineWidth',lw(1))
    hold all
    errorbar(dBins,mean(corrstrucdistOFF(nlv,:,:),3),errOFF(nlv,:),'r','LineWidth',lw(2))
    xlabel('distance bins [mm]')
    ylabel('r_sc')
    legend('Stimulus ON','Stimulus OFF')
    title(['Spatial Structure of r_sc for LV - ' num2str(nlv)])
end
subplot(2,2,4)
errorbar(dBins,mean(rem_corrstrucdistON(1,:,:),3)+0.00001,remerrON,'k','LineWidth',lw(1))
hold all
errorbar(dBins,mean(corrstrucdistON(1,:,:),3)+0.00001,remerrOFF,'c','LineWidth',lw(2))
xlabel('distance bins [mm]')
ylabel('r_sc')
legend('Stimulus ON','Stimulus OFF')
title('Remaining Spatial Structure of r_sc')
suptitle('r_sc along distances. For LVs and Remaining Covs')

% Along angular bins

for nlv = 1:size(corrstrucdistON,1)
    aerrON(nlv,:) = sem(mean(corrstrucdistON(nlv,:,:),2))*ones(1,length(aBins));
    aerrOFF(nlv,:) = sem(mean(corrstrucdistOFF(nlv,:,:),2))*ones(1,length(aBins));
end

aremerrON = sem(mean(rem_corrstrucdistON(1,:,:),2))*ones(1,length(aBins));
aremerrOFF = sem(mean(rem_corrstrucdistOFF(1,:,:),2))*ones(1,length(aBins));

% Plot
figure(2)
for nlv = 1:params.lv
    subplot(2,2,nlv)
    errorbar(rad2deg(aBins),mean(corrstrucdistON(nlv,:,:),2),aerrON(nlv,:),'b','LineWidth',lw(1))
    hold all
    errorbar(rad2deg(aBins),mean(corrstrucdistOFF(nlv,:,:),2),aerrOFF(nlv,:),'r','LineWidth',lw(2))
    xlabel('Angular bins [deg]')
    ylabel('r_sc')
    legend('Stimulus ON','Stimulus OFF')
    title(['Spatial Structure of r_sc for LV - ' num2str(nlv)])
end
subplot(2,2,4)
errorbar(rad2deg(aBins),mean(rem_corrstrucdistON(1,:,:),2),aremerrON,'k','LineWidth',lw(1))
hold all
errorbar(rad2deg(aBins),mean(corrstrucdistON(1,:,:),2),aremerrOFF,'c','LineWidth',lw(2))
xlabel('Angular Bins [deg]')
ylabel('r_sc')
legend('Stimulus ON','Stimulus OFF')
title('Remaining Spatial Structure of r_sc')
suptitle('r_sc along orientations. For LVs and Remaining Covs')

figure(3)

for i = 1:ndatasets
    ax(i)=subplot(2,2,i);
    polar_Image(squeeze(dataset(i).ON.corrstruc(1,:,:)),binCentersON)
    hold all
    polar_Image(squeeze(dataset(i).OFF.corrstruc(1,:,:)),binCentersOFF)
    title('r_sc along distance and angular bins - 1st LV')
end
h=colorbar;
set(h, 'Position', [.8314 .11 .0581 .8150])
for i=1:3
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
end
h.Label.String = 'r_sc';

figure(4)

for i = 1:ndatasets
    ax(i)=subplot(2,2,i);
    polar_Image(squeeze(dataset(i).ON.corrstruc(2,:,:)),binCentersON)
    hold all
    polar_Image(squeeze(dataset(i).OFF.corrstruc(2,:,:)),binCentersOFF)
    title('r_sc along distance and angular bins - 2nd LV')
    
end
h=colorbar;
set(h, 'Position', [.8314 .11 .0581 .8150])
for i=1:3
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
end
h.Label.String = 'r_sc';

figure(5)

for i = 1:ndatasets
    ax(i)=subplot(2,2,i);
    polar_Image(squeeze(dataset(i).ON.rem_corrstruc(1,:,:)),binCentersON)
    hold all
    polar_Image(squeeze(dataset(i).OFF.rem_corrstruc(1,:,:)),binCentersOFF)
    title('r_sc along distance and angular bins - Remaining')
end
h=colorbar;
set(h, 'Position', [.8314 .11 .0581 .8150])
for i=1:3
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
end
h.Label.String = 'r_sc';


figure(6)

ax(1)=subplot(3,1,1);
polar_Image(squeeze(corrstrucdistON(1,:,:)),binCentersON)
hold all
polar_Image(squeeze(corrstrucdistOFF(1,:,:)),binCentersOFF)
title('Avg r_sc along distance and angular bins - 1st LV')

ax(2)=subplot(3,1,2);
polar_Image(squeeze(corrstrucdistON(2,:,:)),binCentersON)
hold all
polar_Image(squeeze(corrstrucdistOFF(2,:,:)),binCentersOFF)
title('Avg r_sc along distance and angular bins - 2nd LV')

ax(3)=subplot(3,1,3);
polar_Image(squeeze(rem_corrstrucdistON(1,:,:)),binCentersON)
hold all
polar_Image(squeeze(rem_corrstrucdistOFF(1,:,:)),binCentersOFF)
title('Avg r_sc along distance and angular bins - Remaining')

h=colorbar;
set(h, 'Position', [.8314 .11 .0581 .8150])
for i=1:3
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
end
h.Label.String = 'r_sc';


%% Variances, weights etc

figure(7)

for i = 1:params.ndatasets
    subplot(2,1,1)
    hold all
    c_on = sort(dataset(i).modelON.C);
    plot([1:size(c_on,1)],c_on(:,1),'b.')
    plot([1:size(c_on,1)],c_on(:,2),'r.')
    xlabel('sorted cells')
    ylabel('weights')
    legend('LV1','LV2')
    title('cells sorted by weights with and without LV - ON')
    subplot(2,1,2)
    hold all
    c_off = sort(dataset(i).modelOFF.C);
    plot([1:size(c_off,1)],c_off(:,1),'b.')
    plot([1:size(c_off,1)],c_off(:,2),'r.')
    xlabel('sorted cells')
    ylabel('weights')
    legend('LV1','LV2')
    title('cells sorted by weights with and without LV - OFF')
end

figure(8)

for i = 1:params.ndatasets
    subplot(2,1,1)
    hold all
    c_on_rem = sort(diag(dataset(i).modelON.rem));
    plot([1:size(c_on_rem,1)],c_on_rem,'.b')
    xlabel('sorted cells')
    ylabel('resid cov')
    legend('Residual Covariances')
    title('cells sorted by residual covariances - ON')
    subplot(2,1,2)
    hold all
    c_off_rem = sort(diag(dataset(i).modelOFF.rem));
    plot([1:size(c_off_rem,1)],c_off_rem,'.r')
    xlabel('sorted cells')
    ylabel('resid cov')
    legend('Residual Covariances')
    title('cells sorted by residual covariances - OFF')
end

figure(9)
for i = 1:params.ndatasets
    subplot(2,1,1)
    hold all
    c_on = (dataset(i).modelON.C);
    plot([1:size(c_on,1)],sort(diag(c_on(:,1)*c_on(:,1)')),'b.')
    plot([1:size(c_on,1)],sort(diag(c_on(:,2)*c_on(:,2)')),'r.')
    xlabel('sorted cells')
    ylabel('cov')
    legend('LV1','LV2')
    title('cells sorted by cov - ON')
    subplot(2,1,2)
    hold all
    c_off = (dataset(i).modelOFF.C);
    plot([1:size(c_off,1)],sort(diag(c_off(:,1)*c_off(:,1)')),'b.')
    plot([1:size(c_off,1)],sort(diag(c_off(:,2)*c_off(:,2)')),'r.')
    xlabel('sorted cells')
    ylabel('cov')
    legend('LV1','LV2')
    title('cells sorted by cov - OFF')
end

figure(10)
for i = 1:params.ndatasets
    subplot(2,1,1)
    hold all
    c_on = sort(dataset(i).modelON.ve);
    plot([1:size(c_on,1)],c_on,'b.')
    xlabel('sorted cells')
    ylabel('variance explained by 2 LV')
    legend('Variance Explained - ON')
    title('cells sorted by VE - ON')
    subplot(2,1,2)
    hold all
    c_off = sort(dataset(i).modelOFF.ve);
    plot([1:size(c_off,1)],c_off,'r.')
    xlabel('sorted cells')
    ylabel('variance explained by 2 LV')
    legend('Variance Explained - OFF')
    title('cells sorted by VE - OFF')
end

%% LVwise hotspots

% Remap weights to the Utah array

for currset = 1:ndatasets
    figure(10+currset)
    
    %SingleColorBar()
    ax(1)=subplot(2,2,1);
    mvlv1on{currset} = image_mapArray(dataset(currset).modelON.C(:,1)', map{currset}, 1);
    xlabel('X-Electrodes')
    ylabel('Y-Electrodes')
    title('1st LV - ON')
    ax(2)=subplot(2,2,2);
    mvlv1off{currset} = image_mapArray(dataset(currset).modelOFF.C(:,1)', map{currset}, 1);
    xlabel('X-Electrodes')
    ylabel('Y-Electrodes')
    title('1st LV - OFF')
    ax(3)=subplot(2,2,3);
    mvlv2on{currset} = image_mapArray(dataset(currset).modelON.C(:,2)', map{currset}, 1);
    xlabel('X-Electrodes')
    ylabel('Y-Electrodes')
    title('2nd LV - ON')
    ax(4)=subplot(2,2,4);
    mvlv2off{currset} = image_mapArray(dataset(currset).modelOFF.C(:,2)', map{currset}, 1);
    xlabel('X-Electrodes')
    ylabel('Y-Electrodes')
    title('2nd LV - OFF')
    
    tit = ['Dataset # : ' num2str(currset)];
    suptitle(tit)
    
    h=colorbar;
    set(h, 'Position', [.8314 .11 .0581 .8150])
    for i=1:4
        pos=get(ax(i), 'Position');
        set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
    end
    h.Label.String = 'Sitewise Correlated Activity';
    
end

%% Timescale box plot?

figure(16)

tausON = [];
tausOFF = [];

for currset = 1:ndatasets
    
    tausON = [tausON;dataset(currset).modelON.tau']; % Tau is in bin size. 50ms is the binsize
    tausOFF = [tausOFF;dataset(currset).modelOFF.tau'];
    
end

tausON = tausON*0.050;
tausOFF = tausOFF*0.050;

tausON = 1./tausON;
tausOFF = 1./tausOFF;

taus = [tausON(:,1) tausOFF(:,1) tausON(:,2) tausOFF(:,2)];

%boxplot(taus,'notch','off','labels',{'LV1 - ON','LV1-OFF','LV2-ON','LV2-OFF'})
ax = gca;
h = notBoxPlot(taus,[1:4],0.1,'sdline');
d = [h.data];
set(d(1:1:end), 'markerfacecolor', [0.4,1,0.4], 'color', [0,0.4,0]);
set(d, 'markersize', 5);
ax.XTickLabel = {'LV1-ON','LV1-OFF','LV2-ON','LV2-OFF'};
ylabel('Timescale in Hz')
xlabel('Latent Variables and Stimulus Conditions')
title('Timescales of the Latent Processes')

%% More plots - Develop this section

% Naive factor loadings vs fitted

figure(17)

for i = 1:ndatasets
    
    Con=(dataset(i).modelON.C);
    nCon = naiveCon{i};
    subplot(2,2,i)
    scatter((Con(:,1)),(nCon(:,1)),2,'b')
    hold all
    scatter((Con(:,2)),(nCon(:,2)),2,'r')
    line(xlim,ylim,'Color','k')
    xlabel('Fitted Factor Loadings for LVs')
    ylabel('Naive Factor Loadings')
    legend('LV 1','LV 2')
    title(['Dataset # : ' num2str(i)])
    tit = ['Naive vs fitted factor loadings'];
    suptitle(tit)
end

figure(18)

for i = 1:ndatasets
    
    Ron=diag(dataset(i).modelON.R);
    nRon = diag(naiveRon{i});
    subplot(2,2,i)
    scatter((Ron(:,1)),(nRon(:,1)),2)
    line(xlim,ylim,'Color','k')
    xlabel('Fitted residCov for LVs')
    ylabel('Private Noise')
    title(['Dataset # : ' num2str(i)])
    tit = ['Noise vs residCov'];
    suptitle('Comparison of noise vs residCov')
end
