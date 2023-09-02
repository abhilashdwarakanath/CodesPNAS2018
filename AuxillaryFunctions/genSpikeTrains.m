function spikeTrains = genSpikeTrains(spikeTimes)

for i = 1:size(spikeTimes,1)
    for j = 1:size(spikeTimes,2)
        st = zeros(10300,1);
        if sum(spikeTimes{i,j})~=0
            st(spikeTimes{i,j}) = 1;
        end
        spikeTrains(:,j,i) = st;
    end
end
