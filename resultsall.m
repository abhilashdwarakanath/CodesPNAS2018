clear all
close all
clc

% This script produces the results and plots them. All four JC CCGs and SC
% CCGs and their patterns.

%% Load data and Analyse

datset = 4;

normfac = [100 30 100 100];

% For stim on

load('SCallON0.mat');
load('JCallON0.mat');

fighands = [1 2 3 4];

for i = 1:datset
    [WIJCon{i},WISCon{i}] = processdatas(resultsJC,resultsSC,i,normfac(i),fighands(i));
end

clear resultsSC
clear resultsJC

% For stim off

load('SCallOff0.mat');
load('JCallOff0.mat');

fighands = [11 12 13 14];

for i = 1:datset
    [WIJCoff{i},WISCoff{i}] = processdatas(resultsJC,resultsSC,i,normfac(i),fighands(i));
end

clear resultsSC
clear resultsJC

%% Get averages of all four sets

% ON

%SC-CCG
for i = 1:2*datset
    for j = 1:datset
        sumccscon{i}(j,:)=WISCon{j}.SC{i};
    end
end

for i=1:2*datset
    ccsconavged(:,i)=mean(sumccscon{i},1);
end

% Rccg

for j = 1:datset % data set
    for k = 1:9 %int wins
        sumrccgon(j,k,:)=WISCon{j}.Rccg(k,:);
    end
end

rccgonavged = mean(sumrccgon,1);
rccgonavged = squeeze(rccgonavged);

% Rccg for different distances

% JC-CCG
for i = 1:2*datset
    for j = 1:datset
        for k = 1:datset
            sumccjcon{i}{k}(j,:)=WIJCon{j}{k}.JC{i};
        end
    end
end

for i=1:2*datset
    for k=1:datset
    ccjconavged(:,i,k)=mean(sumccjcon{i}{k},1);
    end
end

%SC-CCG
for i = 1:2*datset
    for j = 1:datset
        sumccscoff{i}(j,:)=WISCoff{j}.SC{i};
    end
end

for i=1:2*datset
    ccscoffavged(:,i)=mean(sumccscoff{i},1);
end

% Rccg

for j = 1:datset % data set
    for k = 1:9 %int wins
        sumrccgoff(j,k,:)=WISCoff{j}.Rccg(k,:);
    end
end

rccgoffavged = mean(sumrccgoff,1);
rccgoffavged = squeeze(rccgoffavged);

% Rccg for different distances

% JC-CCG
for i = 1:2*datset
    for j = 1:datset
        for k = 1:datset
            sumccjcoff{i}{k}(j,:)=WIJCoff{j}{k}.JC{i};
        end
    end
end

for i=1:2*datset
    for k=1:datset
    ccjcoffavged(:,i,k)=mean(sumccjcoff{i}{k},1);
    end
end

%% Plot averages

ax = linspace(0.25,4.25,8);

pon = polyfit(ax,rccgonavged(7,:),2);
poff = polyfit(ax,rccgoffavged(7,:),2);
fitax = linspace(ax(1),ax(end),length(ax)*10);
fitfunc(1,:) = polyval(pon,fitax);
fitfunc(2,:) = polyval(poff,fitax);

figure(15)
dist = 8;
weight = linspace(0.3,2.4,8);
subplot(421)
x=[1 5 10 20 50 100 500 1000 2500];
for i = 1:dist
    semilogx(x,rccgonavged(:,i),'r-','LineWidth',weight(i))
    hold on
    semilogx(x,rccgoffavged(:,i),'k-','LineWidth',weight(i))
end
xlabel('Integration time windows')
ylabel('Rccg')
title('r_ccgs against Int Wins at different distance bins')
legend('0.25-0.75 On', '0.25-0.75 Off', '0.76-1.25', '0.76-1.25', '1.26-1.75', '1.26-1.75', '1.76-2.25', '1.76-2.25', '2.26-2.75', '2.26-2.75', '2.76-3.25', '2.76-3.25', '3.26-3.75', '3.26-3.75', '3.76-4.25' , '3.76-4.25');

weight = linspace(0.3,2.7,9);
subplot(422)
for i = 1:length(weight)
    plot(ax,rccgonavged(i,:),'r.-','LineWidth',weight(i))
    hold on
    plot(ax,rccgoffavged(i,:),'k.-','LineWidth',weight(i))
end
xlabel('distance bins')
ylabel('Rccg')
title('r_ccgs against distance bins at different Int Wins')
legend('4ms On', '4ms Off', '20ms', '20ms', '40ms', '40ms', '80ms', '80ms', '200ms', '200ms', '400ms', '400ms', '500ms', '500ms', '1000ms', '1000ms', '2000ms', '2000ms');

subplot(423)
plot(ax,rccgonavged(7,:),'-ob');
hold on
plot(ax,rccgoffavged(7,:),'-or','LineWidth',1.5);
plot(fitax,fitfunc(1,:),'.k')
plot(fitax,fitfunc(2,:),'.k','LineWidth',1.5)
xlabel('distance bins')
ylabel('r_ccg')
title('Strength of correlation across distance in the PFC')
legend('Shuffle Corrected - ON','Shuffle Corrected - Off', 'SC-fit - ON','SC-fit - Off')

linespec = {'b','r','g','c','k'};
ax = repmat(ax,5,1);

lagsSC = linspace(-5,5,10001);
lagsJC = linspace(-0.5,0.5,2001);

Smid = (round(length(lagsSC)/2))-1;
Jmid = (round(length(lagsJC)/2))-1;

jitwin=[10 20 50 100];
for i = 1:length(jitwin)
    wind=jitwin(i)/2;
    for k = 1:dist
        AOCon(i,k) = trapz(ccjconavged(Jmid-wind:Jmid+wind,k,i));
    end
end

for k = 1:dist
    AOCon(5,k) = trapz(ccsconavged(Smid-jitwin(3):Smid+jitwin(3),k));
end

for i=1:5
    AOCon(i,:)=AOCon(i,:)/AOCon(i,1);
end

for i = 1:length(jitwin)
    wind=jitwin(i)/2;
    for k = 1:dist
        AOCoff(i,k) = trapz(ccjcoffavged(Jmid-wind:Jmid+wind,k,i));
    end
end

for k = 1:dist
    AOCoff(5,k) = trapz(ccscoffavged(Smid-jitwin(3):Smid+jitwin(3),k));
end

for i=1:5
    AOCoff(i,:)=AOCoff(i,:)/AOCoff(i,1);
end

subplot(424)
for i = 1:5
    plot(ax(i,:),AOCon(i,:),'Color',linespec{i});
    hold on
    plot(ax(i,:),AOCoff(i,:),'Color',linespec{i},'LineWidth',1.5);
end
xlabel('Distance in mm')
ylabel('Normalised AOC under the CCG')
title('Fig 6')
legend('JC-10ms - ON', 'JC-10ms - Off','JC-20ms','JC-20ms','JC-50ms','JC-50ms','JC-100ms','JC-100ms', 'SC', 'SC');

subplot(4,2,[5 6])
for i=1:dist
    
    plot(lagsSC,(ccsconavged(:,i))+(i/75),'k','LineWidth',2)
    hold on
    plot(lagsJC,(ccjconavged(:,i,3))+(i/75),'Color',[0.75 0.75 0.75])
end
xlabel('lags (s)')
ylabel('Coincidence/spk')
title('Shuffle and jitter corrected CCGs compared across distance bins - ON')
legend('Shuffle Corrected CCGs','Jitter Corrected CCGs-50ms')
axis([-5 5 0 0.12])

subplot(4,2,[7 8])
for i=1:dist
    
    plot(lagsSC,(ccscoffavged(:,i))+(i/75),'k','LineWidth',2)
    hold on
    plot(lagsJC,(ccjcoffavged(:,i,3))+(i/75),'Color',[0.75 0.75 0.75])
end
xlabel('lags (s)')
ylabel('Coincidence/spk')
title('Shuffle and jitter corrected CCGs compared across distance bins - Off')
legend('Shuffle Corrected CCGs','Jitter Corrected CCGs-5ms')
axis([-5 5 0 0.12])

figure(151)
for i=1:8
    fullwidthON(i) = fwhm(lagsSC,(ccsconavged(:,i)));
    fullwidthOff(i) = fwhm(lagsSC,(ccscoffavged(:,i)));
    [sON(i), mON(i)]=gaussfit(lagsSC,(ccsconavged(:,i)));
    [sOff(i), mOff(i)]=gaussfit(lagsSC,(ccscoffavged(:,i)));
end

subplot(311)
bar(ax(1,:),[fullwidthON'*1000 fullwidthOff'*1000])
xlabel('Distance bins in mm')
ylabel('Synchrony window in ms')
legend('On','Off')
title('Synchrony windows around 0-lag')

subplot(312)
bar(ax(1,:),[mON'*1000 mOff'*1000])
xlabel('Distance bins in mm')
ylabel('Shift of 0-peak in ms')
legend('On','Off')

subplot(313)
plot(lagsSC,normalise((ccsconavged(:,1))),'-b')
hold on
plot(lagsSC,normalise((ccscoffavged(:,1))),'-r')
mfon = gaussmf(lagsSC,[sON(1) mON(1)]);
mfoff = gaussmf(lagsSC,[sOff(1) mOff(1)]);
plot(lagsSC,mfon,'.k')
plot(lagsSC,mfon,'.k')
xlabel('lags(s)')
ylabel('coincidence/spk')
legend('SC-CCG ON (1000ms Int Win)', 'SC-CCG Off','Gaussian fit - ON','Gaussian fit - Off')
title('SC-CCG data vs fit')

%% Plot together

% weight = linspace(0.3,2.4,8);
% for k = 1:4
% figure(k+20)
% x=[1 5 10 20 50 100 500 1000 2500];
% for i = 1:8
%     semilogx(x,WISCon{k}.Rccg(:,i),'or-','LineWidth',weight(i))
%     hold on
%     semilogx(x,WISCoff{k}.Rccg(:,i),'ok-','LineWidth',weight(i))
% end
% xlabel('Integration time windows')
% ylabel('Rccg')
% title('r_ccgs against Int Wins at different distance bins ON vs OFF - MUA')
% legend('0.25-0.75 On', '0.25-0.75 Off', '0.76-1.25', '0.76-1.25', '1.26-1.75', '1.26-1.75', '1.76-2.25', '1.76-2.25', '2.26-2.75', '2.26-2.75', '2.76-3.25', '2.76-3.25', '3.26-3.75', '3.26-3.75', '3.76-4.25' , '3.76-4.25');
% weight = linspace(0.3,2.7,9);
% figure(k+25)
% for i = 1:size(WISCon{k}.Rccg,1)
%     plot(ax,WISCon{k}.Rccg(i,:),'r.-','LineWidth',weight(i))
%     hold on
%     plot(ax,WISCoff{k}.Rccg(i,:),'k.-','LineWidth',weight(i))
% end
% xlabel('distance bins')
% ylabel('Rccg')
% title('r_ccgs against distance bins at different Int Wins ON vs OFF')
% legend('4ms On', '4ms Off', '20ms', '20ms', '40ms', '40ms', '80ms', '80ms', '200ms', '200ms', '400ms', '400ms', '500ms', '500ms', '1000ms', '1000ms', '2000ms', '2000ms');
% end