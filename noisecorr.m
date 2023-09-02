function [rawncon, rawncoff] = noisecorr(SU,c)

nUnits = size(SU(c).spikesBySampleIndex);

timebinson = 0:500:5000;
timebinsoff = 5001:500:10030;

params.ntrials = 200;
params.fulltime = 5000;

sraon = zeros((length(timebinson)-1),params.ntrials,nUnits(2));
sraoff = zeros((length(timebinson)-1),params.ntrials,nUnits(2));

for i = 1:nUnits(2)
    for j = 1:params.ntrials
        spkson = SU(c).spikesBySampleIndex{i}{j}(SU(c).spikesBySampleIndex{i}{j}<=params.fulltime);
        spksoff = SU(c).spikesBySampleIndex{i}{j}(SU(c).spikesBySampleIndex{i}{j}>params.fulltime);
        
        for b = 1:length(timebinson)-1
            sraon(b,j,i) = numel(spkson(spkson>timebinson(b)&spkson<timebinson(b+1)));
            sraoff(b,j,i) = numel(spksoff(spksoff>timebinsoff(b)&spksoff<timebinsoff(b+1)));
        end
    end
end

zsraon = zeros((length(timebinson)-1),params.ntrials,nUnits(2));
zsraoff = zeros((length(timebinson)-1),params.ntrials,nUnits(2));

czraon = zscore(sraon,1,2);
czraoff = zscore(sraoff,1,2);

for i = 1:length(timebinson)-1
    ccovon(i,:,:) = cov(squeeze(czraon(i,:,:)));
    ccovoff(i,:,:) = cov(squeeze(czraoff(i,:,:)));
end

for i = 1:length(timebinson)-1
    covon(i,:,:) = corrcov(squeeze(ccovon(i,:,:)));
    covoff(i,:,:) = corrcov(squeeze(ccovoff(i,:,:)));
end

% for i = 1:length(timebinson)-1
%     for j = 1:nUnits(2)
%         zsraon(i,:,j) = zscore(sraon(i,:,j));
%         zsraoff(i,:,j) = zscore(sraoff(i,:,j));
%     end
% end
% 
% covon = zeros(length(timebinson)-1,nUnits(2),nUnits(2));
% covoff = zeros(length(timebinson)-1,nUnits(2),nUnits(2));
% 
% for i = 1:nUnits(2)
%     for j = 1:nUnits(2)
%         for k = 1:length(timebinson)-1
%             zsc1on = zsraon(k,:,i);
%             zsc2on = zsraon(k,:,j);
%             
%             zsc1off = zsraoff(k,:,i);
%             zsc2off = zsraoff(k,:,j);
%             
%             covon(k,i,j) = corr(zsc1on',zsc2on');
%             covoff(k,i,j) = corr(zsc1off',zsc2off');
%         end
%     end
% end

%covon = (1/2000)*(zsraon'*zsraon);
%covoff = (1/2000)*(zsraoff'*zsraoff);

rawncon = nanmean(covon,1);
rawncoff = nanmean(covoff,1);

%rawncon = corr(sraon');
%rawncoff = corr(sraoff');

%rawncon = corrcov(covon);
%rawncoff = corrcov(covoff);