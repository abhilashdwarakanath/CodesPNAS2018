function [ncon, ncoff] = shufflednoisecorr(SU,c)

nUnits = size(SU(c).spikesBySampleIndex);

%ncon = zeros(nUnits(2),nUnits(2));
%ncoff = zeros(nUnits(2),nUnits(2));

timebinson = 0:500:5000;
timebinsoff = 5000:500:10030;

for i = 1:nUnits(2)
    for j = 1:nUnits(2)
        
        if i == j
            sraon = zeros(length(timebinson)-1,200);srbon = zeros(length(timebinson)-1,200);
            sraoff = zeros(length(timebinson)-1,200);srboff = zeros(length(timebinson)-1,200);
            for k = 1:200
                spikesa = SU(c).spikesBySampleIndex{i}{k};
                spikesaon = spikesa(spikesa<=5000);
                spikesaoff = spikesa(spikesa>5000);
                spikesb = SU(c).spikesBySampleIndex{j}{k};
                spikesbon = spikesb(spikesb<=5000);
                spikesboff = spikesb(spikesb>5000);
                for b = 1:length(timebinson)-1
                    sraon(b,k) = numel(spikesaon(spikesaon>timebinson(b)&spikesaon<timebinson(b+1)));
                    srbon(b,k) = numel(spikesbon(spikesbon>timebinson(b)&spikesbon<timebinson(b+1)));
                    sraoff(b,k) = numel(spikesaoff(spikesaoff>timebinsoff(b)&spikesaoff<timebinsoff(b+1)));
                    srboff(b,k) = numel(spikesboff(spikesboff>timebinsoff(b)&spikesboff<timebinsoff(b+1)));
                end
            end
            for d = 1:200
                sraon(:,d) = zscore(sraon(:,d)); srbon(:,d) = zscore(srbon(:,d));
                sraoff(:,d) = zscore(sraoff(:,d)); srboff(:,d) = zscore(srboff(:,d));
            end
            corrcon = corrcoef(sraon',srbon');
            corrcoff = corrcoef(sraoff',srboff');
            ncon1 = corrcon(1,2);
            ncoff1 = corrcoff(1,2);
        else
            
            sraon = zeros(length(timebinson)-1,200);srbon = zeros(length(timebinson)-1,200);
            sraoff = zeros(length(timebinson)-1,200);srboff = zeros(length(timebinson)-1,200);
            ncon2 = zeros(1,200); ncoff2 = zeros(1,200);
            for nshuffle = 1:200
                
                shuffle = randperm(200);
                fprintf('Dataset %d. Units %d and %d. Shuffle %d\n',c,i,j,nshuffle)
                for k = 1:200
                    spikesa = SU(c).spikesBySampleIndex{i}{k};
                    spikesaon = spikesa(spikesa<=5000);
                    spikesaoff = spikesa(spikesa>5000);
                    spikesb = SU(c).spikesBySampleIndex{j}{shuffle(k)};
                    spikesbon = spikesb(spikesb<=5000);
                    spikesboff = spikesb(spikesb>5000);
                    for b = 1:length(timebinson)-1
                        sraon(b,k) = numel(spikesaon(spikesaon>timebinson(b)&spikesaon<timebinson(b+1)));
                        srbon(b,k) = numel(spikesbon(spikesbon>timebinson(b)&spikesbon<timebinson(b+1)));
                        sraoff(b,k) = numel(spikesaoff(spikesaoff>timebinsoff(b)&spikesaoff<timebinsoff(b+1)));
                        srboff(b,k) = numel(spikesboff(spikesboff>timebinsoff(b)&spikesboff<timebinsoff(b+1)));
                    end
                end
                for d = 1:200
                    sraon(:,d) = zscore(sraon(:,d)); srbon(:,d) = zscore(srbon(:,d));
                    sraoff(:,d) = zscore(sraoff(:,d)); srboff(:,d) = zscore(srboff(:,d));
                end
                corrcon = corrcoef(sraon',srbon');
                corrcoff = corrcoef(sraoff',srboff');
                ncon2(nshuffle) = corrcon(1,2);
                ncoff2(nshuffle) = corrcoff(1,2);
            end
            ncon3  = mean(ncon2);
            ncoff3  = mean(ncoff2);
        end
        if i == j
            ncon(i,j) = ncon1;
            ncoff(i,j) = ncoff1;
        else
            ncon(i,j) = ncon3;
            ncoff(i,j) = ncoff3;
        end
    end
    
end