
clear all
RootDir = 'path_to_SEF_files_not_aligned';
SubjDir = {'patient1','patient2','...'};

% Paradigm{1} = 'resting_state';
% Paradigm{2} = 'stimulation';
% => flattened levels resting_state vs. stimulation and protocols & parts 

SubjReadyForCleanCheck = 2:9;%[2,5,8,7,6,3,4];

%% ======= 1 ======= get files and clean markers...

for s = 1:length(SubjReadyForCleanCheck)

clear FPOI MergedMrk Header HDeName ICeName hdEEGdur icEEGdur MrkDur Prots

[f,d,fp] = my_recursive_listfiles(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),'');
Prots = unique(d);
Prots = Prots(cellfun(@isempty,regexp(Prots,'electrodes_coordinates'))); % ignore "electrodes_coordinates" folder
Prots = Prots(~strcmp(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),Prots)); % if there are files at the subject's root folder level, ignore paths corresponding to these...

% get file paths
for n = 1:length(Prots)
    [f,d,fp] = my_recursive_listfiles(Prots{n},'.sef');
    if ~isempty(fp)
    FP = cellstr(fp);
    FPOI{n} = FP(cellfun(@isempty,regexp(f,'.(mrk|vrb|xls|xlsx|doc|docx|txt|mat)')));
    %     if size(FPOI{n},1)~=2
    %         warning('Number of .sef file is not 2 in %s\n',Prots{n})
    %     end
    clear ext1 path1 file
    for pne = 1:length(FPOI{n})
        [path1{pne} file1{pne} ext1{pne}] = fileparts(FPOI{n}{pne});
    end
    %======== This checks two things at the same time ========
    % a) there are two files (presumably 1 for icEEG and the other
    % for hdEEG, but this can be checked later)
    % b) file format is .sef for all files (no .eeg, .ep, .eph, ....etc.)
    if any(~strcmp(ext1,'.sef')) || length(FPOI{n})~=2
        error('Number of .sef file is not 2 in %s\n',Prots{n})
    end
    else
        FPOI{n} = {};
    end
end
FPOI = FPOI(~cellfun(@isempty,FPOI));

% Check number of channels to know which is icEEG and which is hdEEG
for n = 1:length(FPOI)
    if exist(strcat(FPOI{n}{1},'.vrb'))==2
        Vrb = read_vrb_file(strcat(FPOI{n}{1},'.vrb'));
        Temp1 = regexp(Vrb(15,:),':','split');
        %         if ~(str2num(strtrim(Temp1{1}{2}))>249)
        NumChan = str2num(strtrim(Temp1{1}{2}));
    else
        [~,hdr] = readsef(FPOI{n}{1},'hdr');
        NumChan = hdr.numchan;
    end
    if NumChan ~=257
        % if it is the case then probably the second file is the hdEEG
        % file, not the first
        TempBKP1 = FPOI{n}{1};
        TempBKP2 = FPOI{n}{2};
        FPOI{n}{1} = TempBKP2;
        FPOI{n}{2} = TempBKP1;
%         warning('\nhdEEG file: %s\nicEEG file: %s',FPOI{n}{1},FPOI{n}{2})

        % if NumChan of first file is always 257, it can be safely considered as
        % the hdEEG file
    end
end

% clean markers
for n = 1:length(FPOI)
    %     [hdEEGmrkBegin, hdEEGmrkEnd, hdEEGmrkLabel] = read_mrk_file(strcat(FPOI{n}{1},'.mrk'));
    %     [icEEGmrkBegin, icEEGmrkEnd, icEEGmrkLabel] = read_mrk_file(strcat(FPOI{n}{2},'.mrk'));
    %     hdEEGmrk = dual_read_mkr(strcat(FPOI{n}{1},'.mrk'));
    %     icEEGmrk = dual_read_mkr(strcat(FPOI{n}{2},'.mrk'));
    %     match_vectors(hdEEGmrkLabel,icEEGmrkLabel,1) % hopefully monotically increasing with steps of 1
    %     match_vectors(hdEEGmrk(:,2),icEEGmrk(:,2),1) % hopefully monotically increasing with steps of 1
    
    %======== Check for various issues ========
    % as documented in "dual_clean_mkr.m" and "Marker_file_check_summary_log.txt"
    MergedMrk{n} = dual_clean_mkr(strcat(FPOI{n}{2},'.mrk'),strcat(FPOI{n}{1},'.mrk'));
    
    %======== Check also sensor labels, which should be consistent, i.e.========
    % a) strictly identical for all hdEEG recordings (at least in principle)
    % b) named coherently in icEEG recordings
    [~,Header{n,1}] = readsef(FPOI{n}{1},'hdr');
    [~,Header{n,2}] = readsef(FPOI{n}{2},'hdr');
    clear Temp2
    for asdf = 1:size(Header{n,1}.sensor_labels,1)
        Temp2 = regexp(Header{n,1}.sensor_labels(asdf,:),'\w+','match');
        HDeName{n}{asdf,1} = Temp2{1};
    end
    clear Temp2
    for asdf = 1:size(Header{n,2}.sensor_labels,1)
        Temp2 = regexp(Header{n,2}.sensor_labels(asdf,:),'\w+','match');
        ICeName{n}{asdf,1} = Temp2{1};
    end
    
    %======== Compare number of time markers and recording duration ========
    % get recording duration by dividing number of datapoints by sampling
    % rate:
    hdEEGdur{n} = Header{n,1}.numsamples/Header{n,1}.srate;
    icEEGdur{n} = Header{n,2}.numsamples/Header{n,2}.srate;
    MrkDur{n} = size(MergedMrk{n},1);
    
end

%% ======= 2 ======= check that everything is fine...
cell2mat(hdEEGdur)
cell2mat(icEEGdur)
cell2mat(MrkDur)

cell2mat(cellfun(@str2double,HDeName,'UniformOutput',0))
diff(cell2mat(cellfun(@str2double,HDeName,'UniformOutput',0)),1,2)
reshape(cellstr(char(cellfun(@char,HDeName,'UniformOutput',0))),257,size(HDeName,2))

SizeTmp = min(cellfun(@length,ICeName));
try
    reshape(cellstr(char(cellfun(@char,ICeName,'UniformOutput',0))),SizeTmp,size(ICeName,2))
catch
    warning('Number of electrodes is not the same from one icEEG file to another!')
end

MergedMrk

%% ======= 3 ======= save quality assessment output for later resampling
save(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)},'QA_cleaning_output.mat'),'FPOI','MergedMrk','Header','HDeName','ICeName','hdEEGdur','icEEGdur','MrkDur')

end

%% ======= 4 ======= check dual_clean_mkr.m logs

[f,d,fp] = my_recursive_listfiles(RootDir,'Marker_file_check_summary_log.mat');
for k = 1:size(fp,1)
    load(fp(k,:))
%     matfile(fp(k,:))
    
    NonNumMarkALL(k,1) = length(NonNumMark{1});
    NonNumMarkALL(k,2) = length(NonNumMark{2});
    MultiOnesALL(k,1) = length(MultiOnes{1});
    MultiOnesALL(k,2) = length(MultiOnes{2});
    DoubletsGapsALL(k,1) = sum(DoubletsGaps{1}==0);
    DoubletsGapsALL(k,2) = sum(DoubletsGaps{2}==0);
    ExtraMarkersDiscardedALL(k,1) = sum(ExtraMarkersDiscarded{1}>1);
    ExtraMarkersDiscardedALL(k,2) = sum(ExtraMarkersDiscarded{2}>1);
    
    icEEGbutNOThdEEGALL(k) = length(icEEGbutNOThdEEG);
    hdEEGbutNOTicEEGALL(k) = length(hdEEGbutNOTicEEG);
    PeriodMismatchALL(k) = length(PeriodMismatch);
    LongPeriodsALL(k) = length(LongPeriods);
    Unwarping(k,1) = sum(ExtraMarkersDiscarded{1}<0);
    Unwarping(k,2) = sum(ExtraMarkersDiscarded{1}<0);
    
%     open(strrep(fp(k,:),'.mat','.txt'))
end
% not very informative, more relevant to look at what remains vs. duration
% of recording, or proportion of kept / rejected markers (periods)...

%% ======= Fixing anti-alias filter alignment issue =======

[~,~,fp2] = my_recursive_listfiles(RootDir,'QA_cleaning_output.mat');

% Anti-alias filter alignment flag (depends on which EGI system was used and whether online or offline correction was applied):
AAFAF = [0,1,0,0,0,0,0,1]; % to correct or not (on new versions of NetStations, anti-alias filtering is corrected online by default, on old versions it had to be applied offline and was done for only some, not all, patients)
Delay4AAFA = [8,36,8,8,8,8,36,36]; % value in milliseconds (depends on system used)

for k = 1:size(fp2,1)
    load(fp2(k,:))
    if AAFAF(k)==0
        for n = 1:size(MergedMrk,2)
            MergedMrk{n}(:,3) = MergedMrk{n}(:,3)+Delay4AAFA(k);
        end
    end
    save(fullfile(fileparts(fp2(k,:)),'QA_cleaning_output_AAFA.mat'),'FPOI','HDeName','Header','ICeName','MergedMrk','MrkDur','hdEEGdur','icEEGdur');
end

%% ======= 6 ======= last checks and resample

[~,~,fp3] = my_recursive_listfiles(RootDir,'QA_cleaning_output_AAFA.mat');
try
    for k = 1:size(fp3,1)
%     for k = 7
        load(fp3(k,:))
        %     matfile(fp(k,:))
        
        filepathALL{k} = FPOI;
        MergedMrkALL{k} = MergedMrk;
        %
        hdEEGdurALL{k} = cell2mat(hdEEGdur);
        icEEGdurALL{k} = cell2mat(icEEGdur);
        MrkDurALL{k} = cell2mat(MrkDur);
        %
        HDeNameALL{k} = HDeName;
        ICeNameALL{k} = ICeName;
        
        fprintf('\n...Doing subject %d out of %d...\n\n',k,size(fp3,1))
        for kk = 1:size(FPOI,2)
            mkr = MergedMrk{kk};
%             mkdir(strrep(fileparts(FPOI{kk}{2}),'1_Raw','2_Resampled'));
%             save(fullfile(strrep(fileparts(FPOI{kk}{2}),'1_Raw','2_Resampled'),'cleaned_merged_markers.mat'),'mkr')
            fprintf('\n...Doing block %d out of %d...\n\n',kk,size(FPOI,2))
%             dual_realign(FPOI{kk}{2},FPOI{kk}{1}, fullfile(strrep(fileparts(FPOI{kk}{2}),'1_Raw','2_Resampled'),'cleaned_merged_markers.mat'), 1000)
        end
        
        %     100./hdEEGdurALL{k}.*MrkDurALL{k}
        %     100./icEEGdurALL{k}.*MrkDurALL{k}
        
        %     open(strrep(fp(k,:),'.mat','.txt'))
    end
    sendolmail('youraddress@example.com','Command successfully terminated','Resampling of dual EEG data OK')
catch ME
    if isa(ME,'MException')
        if ~isempty(ME.stack)
            sendolmail('youraddress@example.com','Error during command execution',['"',ME.message,'" Line ',num2str(ME.stack.line),' in ',ME.stack.file])
            warning(['"',ME.message,'" Line ',num2str(ME.stack.line),' in ',ME.stack.file])
        else
            sendolmail('youraddress@example.com','Error during command execution',ME.message)
            warning(ME.message)
        end
    end
end

%% ======= 7 ======= save patient7 data with same format
clear FPOI MergedMrk Header HDeName ICeName hdEEGdur icEEGdur MrkDur Prots

[f,d,fp] = my_recursive_listfiles(fullfile(RootDir,SubjDir{1}),'');
Prots = unique(d);
Prots = Prots(cellfun(@isempty,regexp(Prots,'electrodes_coordinates'))); % ignore "electrodes_coordinates" folder
Prots = Prots(~strcmp(fullfile(RootDir,SubjDir{1}),Prots)); % if there are files at the subject's root folder level, ignore paths corresponding to these...

for n = 1:length(Prots)
    [f,d,fp] = my_recursive_listfiles(Prots{n},'.sef');
    if ~isempty(fp)
    FP = cellstr(fp);
    FPOI{n} = FP(cellfun(@isempty,regexp(f,'.(mrk|vrb|xls|xlsx|doc|docx|txt|mat)')));
    %     if size(FPOI{n},1)~=2
    %         warning('Number of .sef file is not 2 in %s\n',Prots{n})
    %     end
    clear ext1 path1 file
    for pne = 1:length(FPOI{n})
        [path1{pne} file1{pne} ext1{pne}] = fileparts(FPOI{n}{pne});
    end
    %======== This checks two things at the same time ========
    % a) there are two files (presumably 1 for icEEG and the other
    % for hdEEG, but this can be checked later)
    % b) file format is .sef for all files (no .eeg, .ep, .eph, ....etc.)
    if any(~strcmp(ext1,'.sef')) || length(FPOI{n})~=2
        error('Number of .sef file is not 2 in %s\n',Prots{n})
    end
    else
        FPOI{n} = {};
    end
end
FPOI = FPOI(~cellfun(@isempty,FPOI));

for n = 1:length(FPOI)
    if exist(strcat(FPOI{n}{1},'.vrb'))==2
        Vrb = read_vrb_file(strcat(FPOI{n}{1},'.vrb'));
        Temp1 = regexp(Vrb(15,:),':','split');
        %         if ~(str2num(strtrim(Temp1{1}{2}))>249)
        NumChan = str2num(strtrim(Temp1{1}{2}));
    else
        [~,hdr] = readsef(FPOI{n}{1},'hdr');
        NumChan = hdr.numchan;
    end
    if NumChan ~=257
        % if it is the case then probably the second file is the hdEEG
        % file, not the first
        TempBKP1 = FPOI{n}{1};
        TempBKP2 = FPOI{n}{2};
        FPOI{n}{1} = TempBKP2;
        FPOI{n}{2} = TempBKP1;
        warning('\nhdEEG file: %s\nicEEG file: %s',FPOI{n}{1},FPOI{n}{2})
%         error('%s (supposedly hdEEG file) does not have 257 tracks\n',FPOI{n}{1})
    end
end

fs = 1000;
try
    for n = 1:length(FPOI)
        fprintf('\n...Doing block %d out of %d...\n\n',n,size(FPOI,2))
        
        fprintf('\nReading MicroMed icEEG data...\n\n')
        [EEG des1] = dual_load_sef(FPOI{n}{2});
        labels = strrep(des1.channelnames','aux','');
        Folder = strrep(fullfile(fileparts(FPOI{n}{2}),'icEEG'),'1_Raw','2_Resampled');
        mkdir(Folder);
        save(fullfile(Folder,'icEEG.mat'),'EEG','labels','fs');
        
        fprintf('\nReading EGI hdEEG data...\n\n')
        [EEG des2] = dual_load_sef(FPOI{n}{1});
        labels = des2.channelnames';
        Folder = strrep(fullfile(fileparts(FPOI{n}{1}),'hdEEG'),'1_Raw','2_Resampled');
        mkdir(Folder);
        save(fullfile(Folder,'hdEEG.mat'),'EEG','labels','fs');
    end
    sendolmail('youraddress@example.com','Command successfully terminated','Saving .mat files for patient7 done')
catch ME
    if isa(ME,'MException')
        sendolmail('youraddress@example.com','Error during command execution',ME.message)
    end
end

%% ======= MANUAL FIXES... ======= 
% below are manual fixes for some files that were more difficult to
% automatise than quickly fixing manually errors / inconsistencies...

%% SOLUTION FOR a patient9...

RootFolasdf = 'path_to_raw_data_folder';
% second fix after fix for "zero" time stamps triggers
RootFolasdf = 'path_to_new_raw_data_folder';

[~,~,fp] = my_recursive_listfiles(RootFolasdf,'patientID*.*mrk');
% fp = fp(3:end,:); % only stimulation protocols are problematic
for nasdf = 1:size(fp,1)
    [MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(fp(nasdf,:));
    mkrNEW=nan(length(MarkerLabel),1);
    mkrNEW(:,1)=MarkerTimeStart;
    mkrNEW(:,2)=MarkerTimeEnd;
    for nn = 1:length(MarkerLabel)
        if ~isempty(regexp(MarkerLabel{nn},'gidx'))
            TempExpress = regexprep(regexp(MarkerLabel{nn},'gidx','split'),'"','');
            TempExpress2 = regexp(TempExpress{2},'=','split');
            TempExpress3 = regexp(TempExpress2{2},',','split');
            if ~isempty(TempExpress3{1})
                mkrNEW(nn,3) = str2double(char(TempExpress3{1}));
            end
        else
            mkrNEW(nn,3) = nan; % because it is only the last one in this case,
            % and it has the exact same time stamps as the one before
            % (on top of being the same marker label)
        end
    end
    mkrNEW(isnan(mkrNEW(:,3)),:)=[];
    
    % >>>>> DE-UNWRAP MARKER LABELS !!! <<<<<
    mkrNEW2=nan(size(mkrNEW));
    mkrNEW2(:,1)=mkrNEW(:,1);
    mkrNEW2(:,2)=mkrNEW(:,1);
    % >>>>> DE-UNWRAP MARKER LABELS !!! <<<<<
    for nn = 1:length(mkrNEW(:,1))
        if mod(mkrNEW(nn,3),255)~=0
            mkrNEW2(nn,3)=mod(mkrNEW(nn,3),255);
        else
            mkrNEW2(nn,3)=255; % WORKS BECAUSE THERE IS ONLY 1 PERIOD!!
        end
    end
    
    [Path, Filename, Ext] = fileparts(fp(nasdf,:));
    copyfile(fp(nasdf,:),strcat(Path,filesep,Filename,'_BKP',Ext));
    write_mrk_file_Cartool(fp(nasdf,:),mkrNEW2(:,1),mkrNEW2(:,2),cellstr(num2str(mkrNEW2(:,3))));
end

% for fixing again trigger labels de-synchronization for parts 2 & 4...
ProblematicFiles = cellstr(fp(end-3:end,:));

%% FIXING TRIGGER LABELS DE-SYNCHRONIZATION with patient9 for parts 2, 3 and 4...
% NB: fortunately, we did not stop / pause acquisition with resting-state
% data!...

RootFolasdf = 'path_to_patient9_raw_data';

% ProblematicFiles = {'paths_to_mrk_files'};
[~,~,fp] = my_recursive_listfiles(RootFolasdf,'patient9*.*mrk');
ProblematicFiles = cellstr(fp);
ProblematicFiles = ProblematicFiles(cellfun(@isempty,regexp(ProblematicFiles,'BKP'))); % filter out BKP BKP2 .mrk files, just to be sure... (had to redo many times this step...)

for n = 1:length(ProblematicFiles);
    [MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(ProblematicFiles{n});
    mkrNEW=nan(length(MarkerLabel),1);
    mkrNEW(:,1)=MarkerTimeStart;
    mkrNEW(:,2)=MarkerTimeEnd;
    mkrNEW(:,3)=str2double(MarkerLabel);
    
    if n < 58
        % mkr1
        MagicNumber = 0;
    elseif n > 163
        % part 4's magic number: 171
        MagicNumber = 171;
%     elseif n > 109 && n < 164 %=> because last 2 EGI files from part2 on
%     icEEG are on part3 !!!
    elseif n > 109 && n < 164
        % part 3's magic number: 5
        MagicNumber = 5;
%     elseif  n > 57 && n < 110 %=> because last 2 EGI files from part2 on
%     icEEG are on part3 !!!
    elseif  n > 57 && n < 108
        % part 2's magic number: 49
        MagicNumber = 49;
    elseif n == 108 || n == 109
        MagicNumber = 0; % these two files did not need any correction
    end
    
    % except for two mrk files that do not need any correction...
    
    for nn = 1:length(mkrNEW(:,1))
        if (mkrNEW(nn,3)-MagicNumber)<=0
            mkrNEW(nn,3)=255+(mkrNEW(nn,3)-MagicNumber);
        else
            mkrNEW(nn,3)=(mkrNEW(nn,3)-MagicNumber);
        end
    end
    
%     mkrNEW(isnan(mkrNEW(:,3)),:)=[];
    [Path, Filename, Ext] = fileparts(ProblematicFiles{n});
    copyfile(ProblematicFiles{n},strcat(Path,filesep,Filename,'_BKP2',Ext));
    write_mrk_file_Cartool(ProblematicFiles{n},mkrNEW(:,1),mkrNEW(:,2),cellstr(num2str(mkrNEW(:,3))));
end

%% NB:
% path_to_special_SEF_file
% Last trigger removed (103) because with anti-alias filtering alignment,
% we add 36 ms to triggers and in this case there is less than 36 ms of EEG data after the last trigger.
% And this last trigger is far away from the last stimulations, so it does
% not add anything to keep it... (?)

%% SOLUTION FOR patient6...

ProblematicFile = 'path_to_mrk_for_patient6';
[MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(ProblematicFile);
mkrNEW=nan(length(MarkerLabel),1);
mkrNEW(:,1)=MarkerTimeStart;
mkrNEW(:,2)=MarkerTimeEnd;
for nn = 1:length(MarkerLabel)
    if ~isempty(regexp(MarkerLabel{nn},'DI'))
        TempExpress = regexprep(regexp(MarkerLabel{nn},'DI','split'),'"','');
        TempExpress2 = regexp(TempExpress{2},':','split');
        if strcmp(TempExpress2{1}(1),'N')
            mkrNEW(nn,3) = str2double(char(TempExpress2{1}(2)));
        else
            mkrNEW(nn,3) = str2double(char(TempExpress2{1}));
        end
    else
        mkrNEW(nn,3) = nan; % because it is only the last one in this case,
        % and it has the exact same time stamps as the one before
        % (on top of being the same marker label)
    end
end
mkrNEW(isnan(mkrNEW(:,3)),:)=[];
[Path, Filename, Ext] = fileparts(ProblematicFile);
copyfile(ProblematicFile,strcat(Path,filesep,Filename,'_BKP',Ext));
write_mrk_file_Cartool(ProblematicFile,mkrNEW(:,1),mkrNEW(:,2),cellstr(num2str(mkrNEW(:,3))));

%% SOLUTION FOR other mrk file from patient6

ProblematicFile = 'path_to_other_mrk_file_of_patient6';
[MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(ProblematicFile);
mkrNEW=nan(length(MarkerLabel),1);
mkrNEW(:,1)=MarkerTimeStart;
mkrNEW(:,2)=MarkerTimeEnd;
for nn = 1:length(MarkerLabel)
    if ~isempty(regexp(MarkerLabel{nn},'DI'))
        TempExpress = regexprep(regexp(MarkerLabel{nn},'DI','split'),'"','');
        TempExpress2 = regexp(TempExpress{2},':','split');
        if strcmp(TempExpress2{1}(1),'N')
            mkrNEW(nn,3) = str2double(char(TempExpress2{1}(2)));
        else
            mkrNEW(nn,3) = str2double(char(TempExpress2{1}));
        end
    else
        mkrNEW(nn,3) = nan; % because it is only the last one in this case,
        % and it has the exact same time stamps as the one before
        % (on top of being the same marker label)
    end
end
mkrNEW(isnan(mkrNEW(:,3)),:)=[];
[Path, Filename, Ext] = fileparts(ProblematicFile);
copyfile(ProblematicFile,strcat(Path,filesep,Filename,'_BKP',Ext));
write_mrk_file_Cartool(ProblematicFile,mkrNEW(:,1),mkrNEW(:,2),cellstr(num2str(mkrNEW(:,3))));

%% SOLUTION FOR patient4

RootFolasdf = 'path_to_raw_patient4';

[~,~,fp] = my_recursive_listfiles(RootFolasdf,'patient4*.*mrk');
fp = fp(3:end,:); % only stimulation protocols are problematic
for nasdf = 1:size(fp,1)
    [MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(fp(nasdf,:));
    mkrNEW=nan(length(MarkerLabel),1);
    mkrNEW(:,1)=MarkerTimeStart;
    mkrNEW(:,2)=MarkerTimeEnd;
    for nn = 1:length(MarkerLabel)
        if ~isempty(regexp(MarkerLabel{nn},'DI'))
            TempExpress = regexprep(regexp(MarkerLabel{nn},'DI','split'),'"','');
            TempExpress2 = regexp(TempExpress{2},':','split');
            if ~isempty(TempExpress2{1})
                if strcmp(TempExpress2{1}(1),'N')
                    mkrNEW(nn,3) = str2double(char(TempExpress2{1}(2:end)));
                else
                    mkrNEW(nn,3) = str2double(char(TempExpress2{1}));
                end
            end
        else
            mkrNEW(nn,3) = nan; % because it is only the last one in this case,
            % and it has the exact same time stamps as the one before
            % (on top of being the same marker label)
        end
    end
    mkrNEW(isnan(mkrNEW(:,3)),:)=[];
    [Path, Filename, Ext] = fileparts(fp(nasdf,:));
    copyfile(fp(nasdf,:),strcat(Path,filesep,Filename,'_BKP',Ext));
    write_mrk_file_Cartool(fp(nasdf,:),mkrNEW(:,1),mkrNEW(:,2),cellstr(num2str(mkrNEW(:,3))));
end

%% SOLUTION FOR patient5

ProblematicFile = 'path_to_mrk_patient5';
[MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(ProblematicFile);
mkrNEW=nan(length(MarkerLabel),1);
mkrNEW(:,1)=MarkerTimeStart;
mkrNEW(:,2)=MarkerTimeEnd;
for nn = 1:length(MarkerLabel)
    if ~isempty(regexp(MarkerLabel{nn},'DI'))
        TempExpress = regexprep(regexp(MarkerLabel{nn},'DI','split'),'"','');
        TempExpress2 = regexp(TempExpress{2},':','split');
        if strcmp(TempExpress2{1}(1),'N')
            mkrNEW(nn,3) = str2double(char(TempExpress2{1}(2)));
        else
            mkrNEW(nn,3) = str2double(char(TempExpress2{1}));
        end
    else
        mkrNEW(nn,3) = nan; % because it is only the last one in this case,
        % and it has the exact same time stamps as the one before
        % (on top of being the same marker label)
    end
end
mkrNEW(isnan(mkrNEW(:,3)),:)=[];
[Path, Filename, Ext] = fileparts(ProblematicFile);
copyfile(ProblematicFile,strcat(Path,filesep,Filename,'_BKP',Ext));
write_mrk_file_Cartool(ProblematicFile,mkrNEW(:,1),mkrNEW(:,2),cellstr(num2str(mkrNEW(:,3))));

%% SOLUTION FOR patient1

ProblematicFile = 'path_to_mrk_patient1';
[MarkerTimeStart, MarkerTimeEnd, MarkerLabel] = read_mrk_file(ProblematicFile);
mkrNEW=nan(length(MarkerLabel),1);
mkrNEW(:,1)=MarkerTimeStart;
mkrNEW(:,2)=MarkerTimeEnd;
for nn = 1:length(MarkerLabel)
    if ~isempty(regexp(MarkerLabel{nn},' '))
        TempExpress = regexprep(regexp(MarkerLabel{nn},' ','split'),'"','');
        mkrNEW(nn,3) = str2double(char(TempExpress{2}));
    else
        mkrNEW(nn,3) = nan; % because it is only the last one in this case,
        % and it has the exact same time stamps as the one before
        % (on top of being the same marker label)
    end
end
mkrNEW(isnan(mkrNEW(:,3)),:)=[];
[Path, Filename, Ext] = fileparts(ProblematicFile);
copyfile(ProblematicFile,strcat(Path,filesep,Filename,'_BKP',Ext));
write_mrk_file_Cartool(ProblematicFile,mkrNEW(:,1),mkrNEW(:,2),cellstr(num2str(mkrNEW(:,3))));

%% SOLUTION FOR other mrk of patient1

mkr= dual_read_mkr('path_to_other_mrk_file_of_patient1');
mkrNEW=nan(size(mkr));
mkrNEW(:,1)=mkr(:,1);
mkrNEW(:,2)=mkr(:,1);
% >>>>> DE-UNWRAP MARKER LABELS !!! <<<<<
for nn = 1:length(mkr(:,1))
    if mod(mkr(nn,2),255)~=0
        mkrNEW(nn,3)=mod(mkr(nn,2),255);
    else
        mkrNEW(nn,3)=255; % WORKS BECAUSE THERE IS ONLY 1 PERIOD!! OTHERWISE IT WOULD BE INVALID!
    end
end
copyfile('path_to_other_mrk_file_of_patient1','path_to_other_mrk_file_of_patient1_BKP');
write_mrk_file_Cartool('path_to_other_mrk_file_of_patient1',mkrNEW(:,1),mkrNEW(:,2),cellstr(num2str(mkrNEW(:,3))));


