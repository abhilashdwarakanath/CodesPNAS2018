figure(3)
subplot(2,2,1)
errorbar(dBins,shufflednc.ON(:,1),shufflednc.errON(:,1),'LineWidth',1.5)
hold on
errorbar(dBins,shufflednc.OFF(:,1),shufflednc.errOFF(:,1),'LineWidth',1.5)
errorbar(dBins,nc.ON(:,1),nc.errON(:,1),'LineWidth',2)
errorbar(dBins,nc.OFF(:,1),nc.errOFF(:,1),'LineWidth',2)
xlabel('Distance bins [mm]')
ylabel('noise correlation')
grid on
box off
legend('Shuffled Stim ON','Shuffled Stim OFF','Raw Stim ON','Raw Stim OFF')
title('Makkay 1')
subplot(2,2,2)
errorbar(dBins,shufflednc.ON(:,2),shufflednc.errON(:,2),'LineWidth',1.5)
hold on
errorbar(dBins,shufflednc.OFF(:,2),shufflednc.errOFF(:,2),'LineWidth',1.5)
errorbar(dBins,nc.ON(:,2),nc.errON(:,2),'LineWidth',2)
errorbar(dBins,nc.OFF(:,2),nc.errOFF(:,2),'LineWidth',2)
xlabel('Distance bins [mm]')
ylabel('noise correlation')
grid on
box off
legend('Shuffled Stim ON','Shuffled Stim OFF','Raw Stim ON','Raw Stim OFF')
title('Makkay 2')
subplot(2,2,3)
errorbar(dBins,shufflednc.ON(:,3),shufflednc.errON(:,3),'LineWidth',1.5)
hold on
errorbar(dBins,shufflednc.OFF(:,3),shufflednc.errOFF(:,3),'LineWidth',1.5)
errorbar(dBins,nc.ON(:,3),nc.errON(:,3),'LineWidth',2)
errorbar(dBins,nc.OFF(:,3),nc.errOFF(:,3),'LineWidth',2)
xlabel('Distance bins [mm]')
ylabel('noise correlation')
grid on
box off
legend('Shuffled Stim ON','Shuffled Stim OFF','Raw Stim ON','Raw Stim OFF')
title('Dino 1')
subplot(2,2,4)
errorbar(dBins,shufflednc.ON(:,4),shufflednc.errON(:,4),'LineWidth',1.5)
hold on
errorbar(dBins,shufflednc.OFF(:,4),shufflednc.errOFF(:,4),'LineWidth',1.5)
errorbar(dBins,nc.ON(:,4),nc.errON(:,4),'LineWidth',2)
errorbar(dBins,nc.OFF(:,4),nc.errOFF(:,4),'LineWidth',2)
xlabel('Distance bins [mm]')
ylabel('noise correlation')
grid on
box off
legend('Shuffled Stim ON','Shuffled Stim OFF','Raw Stim ON','Raw Stim OFF')
title('Dino 2')
suptitle('Shuffled and Raw noise corrs along distances - 1s bins - 0.05Contam - Averaging Method')

%% Grand average
figure(4)
errorbar(dBins,nanmean(shufflednc.ON,2),nanmean(shufflednc.errON,2),'LineWidth',2)
hold on
errorbar(dBins,nanmean(shufflednc.OFF,2),nanmean(shufflednc.errOFF,2),'LineWidth',2)
errorbar(dBins,nanmean(nc.ON,2),nanmean(nc.errON,2),'LineWidth',2)
errorbar(dBins,nanmean(nc.OFF,2),nanmean(nc.errOFF,2),'LineWidth',2)
xlabel('Distance bins [mm]')
ylabel('noise correlation')
grid on
box off
legend('Raw Stim ON','Raw Stim OFF','Shuffled Stim ON','Shuffled Stim OFF')
title('Grand Average - 0.05Contam - Averaging Method')