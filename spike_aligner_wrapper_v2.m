% ====== spike_aligner_wrapper_v2 ======

%% Dictionnaries:

% ======== icEEG ========

% ARTEFACTS AND RANDOM WILL NOT BE ALIGNED PERFECTLY, BECAUSE THEY WERE
% ANYHOW JUST ROUGHLY MARKED SO THEY DON'T NEED THAT MUCH TEMPORAL
% PRECISOIN AND 50 ms WON'T MATTER THAT MUCH...

% sub01:
icEEGlist.sub01 = {'"ManSEEG_AG12_neg"'
    '"ManSEEG_AG78_neg"'
    '"ManSEEG_AG78_pos"'
    '"ManSEEG_HAD123"'
    '"ManSEEG_HAD123_neg"'
    '"ManSEEG_HAD12_neg"'
    '"ManSEEG_HAG123_1st"'
    '"ManSEEG_HAG123_2nd"'
    '"ManSEEG_HAG123_neg_pos"'
    '"ManSEEG_HAG123_neg_pos_m"'
    '"ManSEEG_HAG123_pos_neg_m"'
    '"ManSEEG_HAG12_pos"'
    '"ManSEEG_HAG78_neg"'
    '"ManSEEG_HAG78_pos"'
    '"ManSEEG_HPG123"'
    '"ManSEEG_HPG23_neg"'
    '"ManSEEG_HPG78_pos"'
    '"ManSEEG_IPG45_pos"'};

icEEGlist.sub01(:,2) = {'AG1-AG2'
    'AG7-AG8'
    'AG7-AG8'
    'HAD1-HAD2'
    'HAD1-HAD2'
    'HAD1-HAD2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG7-HAG8'
    'HAG7-HAG8'
    'HPG2-HPG3' % before was set to 'HPG1-HPG2', but average clearly shows that HPG2-HPG3 is more relevant
    'HPG2-HPG3'
    'HPG7-HPG8'
    'IPG4-IPG5'};

icEEGlist.sub11 = {'"ManSEEG_0_HAD678"'
    '"ManSEEG_0_HAD_HPD"'
    '"ManSEEG_0_HAG_HPG"'
    '"ManSEEG_0_HPG_pHAG"'
    '"ManSEEG_0_nHAD12"'
    '"ManSEEG_0_pHAD"'
    '"ManSEEG_0_pHAD_pHAG"'
    '"ManSEEG_0_pHAD_pHPD"'
    '"ManSEEG_0_pHAG"'
    '"ManSEEG_0_pHAG_pHPG"'
    '"ManSEEG_0_ppHAD"'
    '"ManSEEG_0_ppHAG"'
    '"ManSEEG_=-pHAG"'};

% 0_HAD678 ("ManSEEG_0_HAD678") is actually not HAD but should rather be HPD678 !!!

icEEGlist.sub11(:,2) = {'HPD5-HPD6'
    'HAD2-HAD3'
    'HAG2-HAG3'
    'HAG2-HAG3'
    'HAD2-HAD3'
    'HAD2-HAD3'
    'HAG2-HAG3'
    'HAD2-HAD3'
    'HAG2-HAG3'
    'HAG2-HAG3'
    'HAD2-HAD3'
    'HAG2-HAG3'
    'HAG2-HAG3'};

icEEGlist.sub12 = {'"ManSEEG_0_HADHPD_l"'
    '"ManSEEG_0_HAD_HPD"'
    '"ManSEEG_0_HAD_HPD_l"'
    '"ManSEEG_0_HAD_HPD_m"'
    '"ManSEEG_0_HPD"'
    '"ManSEEG_0_HPDpHAD"'
    '"ManSEEG_0_TDlat"'
    '"ManSEEG_0_TPSD"'
    '"ManSEEG_0_Tmedlat_simult"'
    '"ManSEEG_0_pTDlat"'};

icEEGlist.sub12(:,2) = {'HPD1-HPD2'
    'HPD1-HPD2'
    'HPD1-HPD2'
    'HPD1-HPD2'
    'HPD1-HPD2'
    'HPD1-HPD2'
    'HPD7-HPD8'
    'TPSD5-TPSD6' % 'TPSD4-TPSD5'
    'HPD1-HPD2' % 'HAD1-HAD2'
    'HPD7-HPD8'};

% '"ManSEEG_0_TDlat"' RELATES TO FIRST POSITIVE, NOT NEGATIVE, DEFLECTION !
% => v3 of spike_aligner now plots the average spike based on the previous
% marking, this should help a lot

icEEGlist.sub33 = {'"ManSEEG_HAD2_HPD2"'
    '"ManSEEG_HAD3_HPD2"'
    '"ManSEEG_HAD3_HPD2m"'
    '"ManSEEG_HAG1_HPG1"'
    '"ManSEEG_HAG1_HPG1_PO"'
    '"ManSEEG_HAG1_HPG1m"'
    '"ManSEEG_HAG1w"'
    '"ManSEEG_HAG2_HPG2"'
    '"ManSEEG_HAG_HPG_TPG"'
    '"ManSEEG_HPD1m"'
    '"ManSEEG_HPD1w"'
    '"ManSEEG_HPG1m"'
    '"ManSEEG_HUP"'};

icEEGlist.sub33(:,2) = {'HAD2-HAD3'
    'HAD2-HAD3'
    'HPD1-HPD2' % 'HAD2-HAD3'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HPD1-HPD2'
    'HPD1-HPD2'
    'HPG1-HPG2'
    'HUP1-HUP2'};

icEEGlist.sub34 = {'"ManSEEG_AG_TPG"'
    '"ManSEEG_AG_TPG_HAG78"'
    '"ManSEEG_HAD123"'
    '"ManSEEG_HAG78"'
    '"ManSEEG_TPG"'
    '"ManSEEG_T_medlat_simult"'
    '"ManSEEG_pAG_TPG"'
    '"ManSEEG_pHAD123"'
    '"ManSEEG_pHAG78"'
    '"ManSEEG_pTPG"'};

icEEGlist.sub34(:,2) = {'ag3-ag4'
    'ag3-ag4'
    'had1-had2'
    'hag7-hag8'
    'tpg1-tpg2'
    'hag7-hag8'
    'ag3-ag4'
    'had1-had2'
    'hag7-hag8'
    'tpg1-tpg2'};

icEEGlist.sub35 = {'"ManSEEG_AD_HAD"'
    '"ManSEEG_HAD"'
    '"ManSEEG_HAD_m"'
    '"ManSEEG_HADm"'
    '"ManSEEG_HADpoly"'
    '"ManSEEG_HAG_HPG"'
    '"ManSEEG_HAG_HPG_pos"'
    '"ManSEEG_HAG_neg"'
    '"ManSEEG_HAG_pos"'
    '"ManSEEG_HAG_pos_m"'
    '"ManSEEG_HPG78_neg"'
    '"ManSEEG_TPDlat"'};
%     '"ManSEEG_HPD_lat"', '"ManSEEG_HPG_lat"' & "ManSEEG_TPG_lat" ignored
% because they consist of a long period with potentially multiple spikes
% but not clear!

icEEGlist.sub35(:,2) = {'AD1-AD2'
    'HAD1-HAD2'
    'HAD1-HAD2'
    'HAD1-HAD2'
    'HAD1-HAD2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HAG1-HAG2'
    'HPG7-HPG8'
    'TPD3-TPD4'
    };

% ======== hdEEG ========

hdEEGlist.sub01 = {'"ManScalp_LT"'};
hdEEGlist.sub01(:,2) = {'T9-TP9'};

hdEEGlist.sub11 = {'"ManScalp_0_t10"'};
hdEEGlist.sub11(:,2) = {'T10-TP10'};

hdEEGlist.sub12 = {'"ManScalp_RT"'}; % '"ManScalp_RT_slow"' & '"ManScalp_TR_slow"' are epochs with duration > 0, so coarse alignment is sufficient
hdEEGlist.sub12(:,2) = {'F8-T8'}; % 'T8-P8' was also a good candidate but the deflection is smoother on the average, so it was likely not used for marking the event and was rather a concomittant deflection with some delay across events with respect to F8-T8

hdEEGlist.sub33 = {''}; % no scalp markers for him!
hdEEGlist.sub33(:,2) = {''};

hdEEGlist.sub34 = {'"ManScalp_LT"'
    '"ManScalp_RT"'};

hdEEGlist.sub34(:,2) = {'T7-P7'
    'T10-TP10'};

hdEEGlist.sub35 = {'"ManScalp_P9"'
    '"ManScalp_T7T9"'
    '"ManScalp_T8T10"'};
hdEEGlist.sub35(:,2) = {'P9-O1'
    'T9-TP9'
    'T10-TP10'};

%% check what hdEEG channels are involved for each spike type (more difficult than icEEG)

FrameLength = 1000;
pat = 2
l = 1
l = 2
f = 1
f = 2
f = 3
OldMarker = regexprep(eval(['hdEEGlist.',Pats{pat},'{l,1}']),'ManScalp_','');
OldFile = [fileparts(regexprep(FileList{f},'2b_Realigned','2_Resampled')),filesep,'hdEEG',filesep,'hdEEG.sef'];
% double banana:
spike_aligner_hdEEG(FileList{f}, eval(['hdEEGlist.',Pats{pat},'{l,1}']), OldFile, OldMarker, eval(['hdEEGlist.',Pats{pat},'{l,2}']), FrameLength ,'double');
% triple banana:
spike_aligner_hdEEG(FileList{f}, eval(['hdEEGlist.',Pats{pat},'{l,1}']), OldFile, OldMarker, eval(['hdEEGlist.',Pats{pat},'{l,2}']), FrameLength ,'triple');


%% Get data and check

% File = spm_select(1,'any','Select .sef EEG files',{},pwd,'.sef$'); % select EEG .sef file
FileList = cellstr(spm_select(Inf,'any','Select .sef EEG files',{},pwd,'.*.sef$')); % select EEG .sef file

Hdr = dual_load_sef_hdr(FileList{1});
[~,~,labelsB] = bipolar_montage(Hdr.channelnames,2);
MarkerOfInterest = {};
for f = 1:length(FileList)
    [Marker_T1,Marker_T2,Marker_Label] = read_mrk_Cartool([FileList{f},'.mrk']);
    UniqueMarkers = unique(Marker_Label);
    if isempty(regexp(FileList{1},'hdEEG', 'once')) && ~isempty(regexp(FileList{1},'icEEG', 'once'))
        MarkerOfInterestTemp = UniqueMarkers(~cellfun(@isempty,regexp(UniqueMarkers,'"ManSEEG_*')));
    elseif ~isempty(regexp(FileList{1},'hdEEG', 'once')) && isempty(regexp(FileList{1},'icEEG', 'once'))
        MarkerOfInterestTemp = UniqueMarkers(~cellfun(@isempty,regexp(UniqueMarkers,'"ManScalp_*')));
    else
        error('Ambiguous filename containing both "icEEG" and "hdEEG"')
    end
    MarkerOfInterest(end+1:end+length(MarkerOfInterestTemp)) = MarkerOfInterestTemp;
    
end; MarkerOfInterest = unique(MarkerOfInterest)';
% Mask = '"ManSEEG_*';
% MaskedMarkers = regexprep(MarkerOfInterest,Mask,'')
labelsB'
% MarkerOfInterest(1:end-1)
MarkerOfInterest


%% Interactive alignment

cprintf('cyan',' /!\\ BEWARE WITH THE POLARITY SWITCH (WAS USUALLY 1st NEGATIVE DEFLECTION) !\n')

FrameLength = 1000;

Pats = fieldnames(icEEGlist); % same for hdEEG, except that for one patient there are no markers on scalp

FileList = cellstr(spm_select(Inf,'any','Select .sef EEG files',{},pwd,'.*.sef$')); % select EEG .sef file

% %% TODO:
% f = 1;
% OldMarker = regexprep(eval(['icEEGlist.',Pats{pat},'{l,1}']),'ManSEEG_','');
% OldFile = [fileparts(regexprep(FileList{f},'2b_Realigned','2_Resampled')),filesep,'icEEG',filesep,'icEEG.sef'];
% [Shifts, ToCheck] = spike_aligner_v3(FileList{f}, eval(['icEEGlist.',Pats{pat},'{l,1}']), OldFile, OldMarker, eval(['icEEGlist.',Pats{pat},'{l,2}']), FrameLength );

Root = 'E:\dualEEG\2b_Realigned';
PatFol = {'patient1'
    'patient11'
    'patient12'
    'patient33'
    'patient34'
    'patient35'};

%% interactive review for icEEG
for pat = 1:length(Pats)
    
    [~,~,FileList] = my_recursive_listfiles(fullfile(Root,PatFol{pat},'resting_state'),'^icEEG.*.sef$');
    FileList = cellstr(FileList);

    for f = 1:length(FileList)
        % backup all .mrk files at once
        copyfile([FileList{f},'.mrk'],spm_file([FileList{f},'.mrk'],'suffix','_BKP_before_align'));
    end
    
%     for l = 15:18 % 8:18 %1:size(eval(['icEEGlist.',Pats{pat}]),1) % that was for sub01
%     for l = 7:13 % 2:size(eval(['icEEGlist.',Pats{pat}]),1) % that was for sub11
%     for l = 9:10 % 8:10 % 1:size(eval(['icEEGlist.',Pats{pat}]),1) % that was for sub12
%     for l = 4:13 % 3:13 % 1:size(eval(['icEEGlist.',Pats{pat}]),1) % that was for sub33
%     for l = 8:10 % 4:10 % 1:size(eval(['icEEGlist.',Pats{pat}]),1) % that was for sub34
%     for l = 2:12 % 1:size(eval(['icEEGlist.',Pats{pat}]),1) % that was for sub34
    for l = 1:size(eval(['icEEGlist.',Pats{pat}]),1)
        
        fprintf('Spike %d/%d...\n',l,size(eval(['icEEGlist.',Pats{pat}]),1))
        
        if ~isempty(eval(['icEEGlist.',Pats{pat},'{l,1}']))
            
            for f = 1:length(FileList)
                
                
                if strcmp(eval(['icEEGlist.',Pats{pat},'{l,1}']),'"ManSEEG_0_TDlat"')
                    cprintf('cyan',' /!\\ FOR THIS SPIKE, IT WAS THE FIRST POSITIVE, NOT NEGATIVE DEFLECTION !\n')
                end
                
                OldMarker = regexprep(eval(['icEEGlist.',Pats{pat},'{l,1}']),'ManSEEG_','');
                OldFile = [fileparts(regexprep(FileList{f},'2b_Realigned','2_Resampled')),filesep,'icEEG',filesep,'icEEG.sef'];
                [Shifts, ToCheck] = spike_aligner_icEEG(FileList{f}, eval(['icEEGlist.',Pats{pat},'{l,1}']), OldFile, OldMarker, eval(['icEEGlist.',Pats{pat},'{l,2}']), FrameLength );
                
                [T1,T2,L] = read_mrk_Cartool([FileList{f},'.mrk']);
                T1(strcmp(L,eval(['icEEGlist.',Pats{pat},'{l,1}']))) = T1(strcmp(L,eval(['icEEGlist.',Pats{pat},'{l,1}'])))-Shifts;
                T2(strcmp(L,eval(['icEEGlist.',Pats{pat},'{l,1}']))) = T2(strcmp(L,eval(['icEEGlist.',Pats{pat},'{l,1}'])))-Shifts;
                write_mrk_file_Cartool([FileList{f},'.mrk'],T1,T2,L);
                
                icEEGshiftsALL{pat,l,f} = Shifts; %#ok<SAGROW>
                icEEG2CheckALL{pat,l,f} = ToCheck; %#ok<SAGROW>
            end
            
        else
            warning('No marker for patient %s!',Pats{pat})
        end
    end
    
end

%% interactive review for hdEEG
for pat = 1:length(Pats)
    
    [~,~,FileList] = my_recursive_listfiles(fullfile(Root,PatFol{pat},'resting_state'),'^hdEEG.*.sef$');
    FileList = cellstr(FileList);
    
    for f = 1:length(FileList)
        % backup all .mrk files at once
        copyfile([FileList{f},'.mrk'],spm_file([FileList{f},'.mrk'],'suffix','_BKP_before_align'));
    end
    
    for l = 1:size(eval(['hdEEGlist.',Pats{pat}]),1)
        
        fprintf('Spike %d/%d...\n',l,size(eval(['hdEEGlist.',Pats{pat}]),1))
        
        if ~isempty(eval(['hdEEGlist.',Pats{pat},'{l,1}']))
            
            for f = 1:length(FileList)
                
                OldMarker = regexprep(eval(['hdEEGlist.',Pats{pat},'{l,1}']),'ManScalp_','');
                OldFile = [fileparts(regexprep(FileList{f},'2b_Realigned','2_Resampled')),filesep,'hdEEG',filesep,'hdEEG.sef'];
                % double banana:
                [Shifts, ToCheck] = spike_aligner_hdEEG(FileList{f}, eval(['hdEEGlist.',Pats{pat},'{l,1}']), OldFile, OldMarker, eval(['hdEEGlist.',Pats{pat},'{l,2}']), FrameLength ,'double');
                % % triple banana:
                % spike_aligner_hdEEG(FileList{f}, eval(['hdEEGlist.',Pats{pat},'{l,1}']), OldFile, OldMarker, eval(['hdEEGlist.',Pats{pat},'{l,2}']), FrameLength ,'triple');
                
                [T1,T2,L] = read_mrk_Cartool([FileList{f},'.mrk']);
                T1(strcmp(L,eval(['hdEEGlist.',Pats{pat},'{l,1}']))) = T1(strcmp(L,eval(['hdEEGlist.',Pats{pat},'{l,1}'])))-Shifts;
                T2(strcmp(L,eval(['hdEEGlist.',Pats{pat},'{l,1}']))) = T2(strcmp(L,eval(['hdEEGlist.',Pats{pat},'{l,1}'])))-Shifts;
                write_mrk_file_Cartool([FileList{f},'.mrk'],T1,T2,L);
                
                hdEEGshiftsALL{pat,l,f} = Shifts; %#ok<SAGROW>
                hdEEG2CheckALL{pat,l,f} = ToCheck; %#ok<SAGROW>
            end
            
        else
            warning('No marker for patient %s!',Pats{pat})
        end
    end
    
end

% DO NOT rename markers with channel names (bipolar), because otherwise we
% loose the amplitude information for example, like in pHAD vs HAD!

%% FIX WRONG CODE ABOVE!!!
Pats{7} = 'sub32';
Pats{8} = 'sub36';
PatFol{7} = 'patient32';
PatFol{8} = 'patient36';
for pat = 1:length(Pats)
    
    [~,~,FileList] = my_recursive_listfiles(fullfile(Root,PatFol{pat},'resting_state'),'^(hd|ic)EEG.*.sef$');
    FileList = cellstr(FileList);
    
    for f = 1:length(FileList)        
        [T1,T2,L] = read_mrk_Cartool(spm_file([FileList{f},'.mrk'],'suffix','_BKP_OK'));
        T1 = T1(cellfun(@isempty,regexp(L,'(EprimeMicromed_|EprimeNetStation_)')));
        T2 = T2(cellfun(@isempty,regexp(L,'(EprimeMicromed_|EprimeNetStation_)')));
        L = L(cellfun(@isempty,regexp(L,'(EprimeMicromed_|EprimeNetStation_)'))); % THIS WAS WRONGLY PLACED BEFORE!
        L = regexprep(L,'(ManSEEG_|ManScalp_)','');
        if ~isempty(L)
            write_mrk_file_Cartool([FileList{f},'.mrk'],T1,T2,L);
        else
            warning('No markers left for %s !',FileList{f})
        end
        
    end
    
end

%% debug

f = 1
f = 2

[Marker_T1,Marker_T2,Marker_Label] = read_mrk_Cartool([FileList{f},'.mrk']);

ThisMarker = '"ManSEEG_HAD2_HPD2"';
ThisMmatch = match_vectors({ThisMarker},Marker_Label,1);
if iscell(ThisMmatch) && cellfun(@isempty,ThisMmatch)
    t = Inf;
    warning('Choose another file, could not locate marker!')
elseif iscell(ThisMmatch)
    t = ThisMmatch{1}(1);
elseif isnumeric(ThisMatch)
    t = ThisMmatch(1);
else
    error('???')
end

clipboard('copy',num2str(Marker_T1(t)))
system(['"C:\Program Files (x86)\Cartool\Cartool.exe" "',FileList{f},'" &'])
% Then press CTRL + G in Cartool and paste clipboard content to go to
% marker !

t = ThisMmatch{1}(2);clipboard('copy',num2str(Marker_T1(t)))
t = ThisMmatch{1}(3);clipboard('copy',num2str(Marker_T1(t)))
