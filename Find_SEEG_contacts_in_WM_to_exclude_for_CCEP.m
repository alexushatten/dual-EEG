%% find SEEG contacts not in GM to exclude for CCEP:

WMseeg = importdata('E:\dualEEG\patient50\T1.hdr.mrk',' ');

[x,y,z,name,ClusterName,Nclus,Fullname] = read_els_file('E:\dualEEG\patient50\de171.els');
ALLseeg = [x,y,z];

IdxWM = match_vectors(cellstr(num2str(WMseeg)),cellstr(num2str(ALLseeg)),1);

% to exclude:
sort(Fullname(IdxWM)')
Fullname(dnif(IdxWM,size(ALLseeg,1)))'

% to keep:
Fullname(~dnif(IdxWM,size(ALLseeg,1)))'

%% now to make pairs to stimulate, we need to account for SEEG contacts in
% WM but next to SEEG contacts in GM, so we apply a dilatation filter:

GMstrict = ~dnif(IdxWM,size(ALLseeg,1));
GMloose = conv(GMstrict+0,[1 1 1])>0;
GMloose =GMloose(2:end-1);

figure; imagesc([GMstrict,GMloose])

[Bipoles,ElecMatch,labelsB,labelsB2] = bipolar_montage(Fullname(GMloose)',1);

% number of minutes of stimulation with 20 repetitions and 5 seconds in
% between each block
size(labelsB2,1)*(20+5)/60
size(labelsB2,1)*(30+5)/60 % 30 repetitions
size(labelsB2,1)*(25+5)/60 % 25 repetitions

