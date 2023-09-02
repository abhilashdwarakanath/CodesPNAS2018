function utahMaps = generatUtahMaps(data)
% utahMaps = generatUtahMaps
% this function generatin Utah maps ONLY for SPIKING channels <# spiking channels = 93>
% 
% <# spiking channels = 93>
% MAINLY FOR 93 CHANNELS with SPIKING activity 
% generating a matrix with 3 column; each column contain on type of indices
% coloums:
% board#    LFP#    spk#

%% MUA
% and partially LFPs
nMUA = numel(data.exp.tetids);

channelsTabel = nan(nMUA, 3);
for iCh = 1 : nMUA
    % pick up usable channels. 'data.exp.tetids' contain usable channels
    tmpIndex = find(data.lfp.params.ttIds == data.exp.tetids(iCh));
    
    channelsTabel(iCh, 1) = data.exp.tetids(iCh);
    channelsTabel(iCh, 2) = tmpIndex;
    channelsTabel(iCh, 3) = iCh;   
end

% Board map (type I)
% which tell each channel number (out of 120) located in which location of real ellectrode
% array
% i. e. anatomically you know each channel will mapped to postion of
% electord, which other channels are it neighbours
% this matrix will be used for function "MUAimage"
% boardMap = ... 
%     [17       20      75      76      101     113     107     108     112     118; ...
% 19      22      77      78      88      111     109     119     116     120; ...
% 21      24      79      80      82      2       115     65      67      66; ...
% 23      34      81      83      84      86      6       117     69      68; ...
% 33      36      85      87      3       8       4       14      71      70; ...
% 35      NaN     42      1       5       7       10      NaN     NaN     NaN; ...
% 37      40      44      46      47      49      9       11      13      NaN; ...
% 39      41      43      45      48      50      51      54      NaN     NaN; ...
% 98      NaN     99      100     102     104     52      53      NaN     55; ...
% 18      73      74      103     105     NaN     106     110     114     NaN];

% LFP map (type II)
% similar to type I, the ony differnece is instead of indicating the #board
% indicated based on 96 recorded channels
% 96 lfp MUA Map - all channels which lfp recorded from them 
% (there is 3 channel there which didn't have spiking actvity)
LFPmap = ...
    [17       20      59      60      77      89      83      84      88      94
19      22      61      62      72      87      85      95      92      96
21      24      63      64      66      2       91      49      51      50
23      26      65      67      68      70      6       93      53      52
25      28      69      71      3       8       4       14      55      54
27      30      34      1       5       7       10      12      16      56
29      32      36      38      39      41      9       11      13      NaN
31      33      35      37      40      42      43      46      15      NaN
74      73      75      76      78      80      44      45      48      47
18      57      58      79      81      NaN     82      86      90      NaN];

% MUA map (type III)
% similar to type I and II, the ony differnece is instead of indicating the
% #board of #lfp, here it's indicated based on 93 channels in which spiking
% activity was recorded
% so 93 MUA Map - only channels contain spiking activity
MUAmap = nan(size(LFPmap));
[MUArow MUAcol] = channelLocFinder((1 : nMUA), channelsTabel, LFPmap);
idx = sub2ind(size(MUAmap), MUArow, MUAcol);
MUAmap(idx) = channelsTabel(:, 3);

% producing a channel Map for later usage - ONLY FOR 93 CHANNELS with SPIKING activity 
% this integrate all I did in "Channel Mapping" section.
% coloums:
% board#    LFP#    spk#    MUArow#     MUAcolumn#

MUAintegratedTabel = nan(nMUA, 5);


MUAintegratedTabel (: , 1 : 3) = channelsTabel;
MUAintegratedTabel(:, 4) = MUArow';
MUAintegratedTabel(:, 5) = MUAcol';

%% SUA
% producing a channel tabel for later usage - ONLY FOR 126 SUA with ISOLATED SPIKING activity 
% coloums:
% board#    LFP#    spk#    MUArow#     MUAcolumn# contaminationValue
nSUA = numel(data.susTT);
SUAintegratedTabel(:, 1) = data.susTT;
SUAintegratedTabel(:, 6) = data.contam(:, 1);

for iUnit = 1 : nSUA
    [tmp_r ~]= find(MUAintegratedTabel(:, 1) == SUAintegratedTabel(iUnit, 1));
    SUAintegratedTabel(iUnit, [2 3 4 5]) = ...
        MUAintegratedTabel(tmp_r, [2 3 4 5]);
end

% producing the SUA map (similar to what we had for LFPmap nad MUAmao)
% ~allocating memory
[nRow_array, nCol_array] = size(LFPmap);
SUAmap{nRow_array, nCol_array} = [];

% channels in which the SUA was detected
% (by choosing first column we are picking the board#)
refColumn = 1;
% ** you can easily use other coumun rather than board# i.e. refColumn = 1 
uniqueSUA = unique(SUAintegratedTabel(:, refColumn)); 

for iuChSUA = [uniqueSUA'] % iuChSUA -> i uniqe channel SUA
    %
    [tmp_index_SUA, ~] = find(SUAintegratedTabel(:, refColumn) == iuChSUA);
    % find the MUArow which correspoond to channel from which SUAs with idex tmp_r_SUA are recorded 
    [tmp_index_MUA, ~] = find(MUAintegratedTabel(:, refColumn) == iuChSUA);
    % pick the location of 'iuChSUA' from MUAintegratedMap
%     [tmp_r_MUA tmp_c_MUA] = MUAintegratedTabel(tmp_index_MUA, [4 5]);
    tmp_loc_MUA = MUAintegratedTabel(tmp_index_MUA, [4 5]);
    SUAmap{tmp_loc_MUA(1), tmp_loc_MUA(2)} = tmp_index_SUA;
end

%% a nice output
% ONE structure containing all the maps
utahMaps.LFP.map = LFPmap;
utahMaps.MUA.channelsTabel = channelsTabel;
utahMaps.MUA.integratedTabel = MUAintegratedTabel;
utahMaps.MUA.map = MUAmap;
utahMaps.SUA.integratedTabel = SUAintegratedTabel;
utahMaps.SUA.map = SUAmap;

