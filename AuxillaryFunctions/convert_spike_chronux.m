function spikeTrains_chronuxFormat = convert_spike_chronux(spikeTrains)
% spikeTrains_chronuxFormat = convert_spike_chronux(spikeTrains)
% author: Natalie Zhavoronkov
%  start: 2013/03/15
% 
% IN:
%   spikeTrains : e.g. data.muSpikesByTime or data.spikesByTime
% OUT:
%   spikeTrains_chronuxFormat : structure -> spikeTrains_chronuxFormat(x).condition(x).trial(x).times

% nSpikes = 
nTrials = size(spikeTrains, 2);
nChannels = size(spikeTrains, 1);

spikeTrains_chronuxFormat = struct();
for iCh = 1 : nChannels
%     for iCond = 1:size(spikeTrains,2)
%         nrTrials = size(spikeTrains{iSpk,iCond},2);
        for iTrial = 1 : nTrials
%             tmp = spikeTrains{iCh, iTrial};
            spikeTrains_chronuxFormat(iCh).trial(iTrial).times = ...
                spikeTrains{iCh, iTrial}(:, 1);
        end
%     end
end

end
