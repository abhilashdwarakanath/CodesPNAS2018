clear all
clc
close all

% This is a program to analyse the coherence data

%% For SUA

phase = 1; % 1=On or 0=off

if phase == 2
    load MUallONSSC310150.mat
else
    load MUallOffSSC3p01.mat
end

f = linspace(10,150,2294); % Generate a list of frequencies between [10 55]Hz.
ker = gausswin(15);
% Average across datasets and distance bins
for i = 1:8
    meanssc(:,i)=(mean(resultsSSC{1}{i},2)+mean(resultsSSC{2}{i},2)+mean(resultsSSC{3}{i},2)+mean(resultsSSC{4}{i},2))/4;
    meanssc(:,i)=conv(meanssc(:,i),ker,'same');
end

% Distance bins (0.5mm)
ax = linspace(0.25,4.25,8);

% Do analysis for freq bins (3 Hz)

f_vals = 10:30:150; %value to find (frequencies in Hz)

stepsize = round(length(f)/max(f));

for i = 1:length(f_vals)-1
    
    tmp1 = abs(f-f_vals(i));
    tmp2 = abs(f-f_vals(i+1));
    
    [idx1 idx1] = min(tmp1); %index of closest value
    [idx2 idx2] = min(tmp2); %index of closest value of the end of the bin
    
    f_int1 = f(idx1); %closest value
    f_int2 = f(idx2); %closest value
    
    % compute for distance bins
    
    for j = 1:length(ax)
        
        coh_est_dist(i,j)= mean(meanssc(idx1:idx2,j));
        
    end
    
end

%% Plot
figure(1)
for i = 2:4
    subplot(3,2,i-1)
    plot(ax,coh_est_dist(i-1,:),'-ok')
    xlabel('distance bins [mm]')
    ylabel('Coherence Estimate')
    title(['Coherence in frequency band : ' num2str(f_vals(i-1)) '-' num2str(f_vals(i)) 'Hz'])
end

subplot(3,2,4)
plot(ax,coh_est_dist(4,:),'-ok')
xlabel('distance bins [mm]')
ylabel('Coherence Estimate')
title(['Coherence in frequency band : ' num2str(100) '-' num2str(130) 'Hz'])

subplot(3,2,[5 6])
mcohestdist = mean(coh_est_dist,1);
plot(ax,mcohestdist,'-or','LineWidth',2)
xlabel('distance bins [mm]')
ylabel('Coherence Estimate')
title('Averaged Coherence Estimate across [10,150] Hz')


figure(2)
imagesc(f,ax,meanssc')
title('Coherence across distances across all frequencies')
xlabel('Frequencies [Hz]')
ylabel('Distances [mm]')