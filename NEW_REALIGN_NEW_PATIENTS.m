% NEW REALIGN, 2018-12-05

%% Downsample icEEG (not hdEEG)

clear all
RootDir = 'E:\dualEEG\1_Raw';
SubjDir = {'patient50','patient48',...
    'patient82','patient96'};

% Anti-alias filter alignment flag:
AAFAF = [0,0,0,0];
Delay4AAFA = [36,36,36,36];

% Paradigm{1} = 'resting_state';
% Paradigm{2} = 'stimulation';
% => flattened levels resting_state vs. stimulation and protocols & parts

SubjReadyForCleanCheck = 1:4; % [1,3,4]; % 1:4; % 3:4;

hd_EEG_sampling_frequency = 1000;

for s = 1:length(SubjReadyForCleanCheck)
    
    clear FPOI MergedMrk Header HDeName ICeName hdEEGdur icEEGdur MrkDur Prots
    
    [f,d,fp] = my_recursive_listfiles(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),'');
    Prots = unique(d);
    Prots = Prots(cellfun(@isempty,regexp(Prots,'electrodes_coordinates'))); % ignore "electrodes_coordinates" folder
    Prots = Prots(~strcmp(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),Prots)); % if there are files at the subject's root folder level, ignore paths corresponding to these...
    
    Prots = Prots(cellfun(@isempty,regexp(Prots,'DO_NOT_USE'))); % IGNORE ALSO FOLDERS MATCHING THIS PATTERN
    
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
    
    for n = 1:length(Prots)
        
        % ======= downsample icEEG traces =======
        %         [hd_EEG,hd_Hdr] = dual_load_sef(FPOI{n}{1});
        [ic_EEG,ic_Hdr] = dual_load_sef(FPOI{n}{2});
        
        %         hd_Nchan = size(hd_EEG,1);
        %         hd_Ntf = size(hd_EEG,2);
        ic_Nchan = size(ic_EEG,1);
        ic_Ntf = size(ic_EEG,2);
        
        ic_EEGi = zeros(ic_Nchan,floor(ic_Ntf/(ic_Hdr.samplingfreq/hd_EEG_sampling_frequency)));
        for c = 1:ic_Nchan
            ic_EEGi(c,:) = interp1(1:ic_Ntf,ic_EEG(c,:),linspace(1,ic_Ntf,floor(ic_Ntf/(ic_Hdr.samplingfreq/hd_EEG_sampling_frequency))));
        end
        
        % ======= downsample triggers timings =======
        try
            [ic_T1,ic_T2,ic_L] = read_mrk_file([FPOI{n}{2},'.mrk']);
        catch
            try
                [ic_T1,ic_T2,ic_L] = read_mrk_Cartool([FPOI{n}{2},'.mrk']);
            catch
                [ic_T1,ic_T2,ic_L] = read_mrk_Excel_export([FPOI{n}{2},'.mrk']);
            end
        end
        ic_T1_new = round(ic_T1/(ic_Hdr.samplingfreq/hd_EEG_sampling_frequency));
        ic_T2_new = round(ic_T2/(ic_Hdr.samplingfreq/hd_EEG_sampling_frequency));
        
        % ======= write output =======
        TargetDir = fileparts(strrep(FPOI{n}{2},'1_Raw','2a_Downsampled'));
        mkdir(TargetDir)
        
        save(fullfile(TargetDir,'icEEGi.mat'),'ic_EEGi','ic_Hdr');
        save(fullfile(TargetDir,'icEEGi_mrk.mat'),'ic_T1_new','ic_T2_new','ic_L','ic_T1','ic_T2');
        
    end
    
end

%% Align icEEG and hdEEG, adjust and relabel markers

for s = 1:length(SubjReadyForCleanCheck)
    
    fprintf('Doing subject %s/%s...\n',num2str(s),num2str(length(SubjReadyForCleanCheck)))
    clear FPOI MergedMrk Header HDeName ICeName hdEEGdur icEEGdur MrkDur Prots
    
    [f,d,fp] = my_recursive_listfiles(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),'');
    Prots = unique(d);
    Prots = Prots(cellfun(@isempty,regexp(Prots,'electrodes_coordinates'))); % ignore "electrodes_coordinates" folder
    Prots = Prots(~strcmp(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),Prots)); % if there are files at the subject's root folder level, ignore paths corresponding to these...
    
    Prots = Prots(cellfun(@isempty,regexp(Prots,'DO_NOT_USE'))); % IGNORE ALSO FOLDERS MATCHING THIS PATTERN
    
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
    
    % Select simultaneous times and resample marker timings
    for n = 1:length(Prots)
        
        fprintf('Doing part %s/%s...\n',num2str(n),num2str(length(Prots)))
        
        % ======= load icEEG traces =======
        fprintf('Loading hdEEG data...\n')
        [hd_EEG,hd_Hdr] = dual_load_sef(FPOI{n}{1});
        
        TargetDir = fileparts(strrep(FPOI{n}{2},'1_Raw','2a_Downsampled'));
        load(fullfile(TargetDir,'icEEGi.mat'),'ic_EEGi','ic_Hdr')
        
        % % Sanity check but systems clocks might be slightly shifted...
        
        % Hdr_ic.day
        % Hdr_ic.hour
        % Hdr_ic.minute
        % Hdr_ic.second
        % Hdr_ic.msecond
        
        % Hdr_hd.day
        % Hdr_hd.hour
        % Hdr_hd.minute
        % Hdr_hd.second
        % Hdr_hd.msecond
        
        % ===========================================
        % ======= downsample triggers timings =======
        % ===========================================
        
        [ic_T1,ic_T2,ic_L] = read_mrk_Cartool([FPOI{n}{2},'.mrk']);
        [hd_T1,hd_T2,hd_L] = read_mrk_Cartool([FPOI{n}{1},'.mrk']);
        
        % % Check that duration is 0 !
        if ~all((ic_T1-ic_T2)==0)
            cprintf('green_','Some icEEG markers have duration > 0 !\n')
            Temp = unique(ic_L((ic_T1-ic_T2)~=0));
            for asdf = 1:length(Temp)
                cprintf('cyan',[Temp{asdf},'\n'])
            end
        end
        if ~all((hd_T1-hd_T2)==0)
            cprintf('green_','Some icEEG markers have duration > 0 !\n')
            Temp = unique(hd_L((hd_T1-hd_T2)~=0));
            for asdf = 1:length(Temp)
                cprintf('cyan',[Temp{asdf},'\n'])
            end
        end
        
        % % Get proper 1st marker:
        warning off %#ok<WNOFF>
        CleanedMrk = dual_clean_mrk_NEW([FPOI{n}{2},'.mrk'],[FPOI{n}{1},'.mrk'],true,255,DINtype);
        warning on %#ok<WNON>
        
        if ~isempty(CleanedMrk) % sometimes there are not overlapping markers !
            
            % ==================================================
            % ======= Anti-Alias Filter Alignment (AAFA) =======
            % ==================================================
            
            DelayAAFA = AAFAF(SubjReadyForCleanCheck(s))*Delay4AAFA(SubjReadyForCleanCheck(s));
            % marker timings hdEEG (CleanedMrk(:,3)) = marker timings hdEEG (CleanedMrk(:,3)) + DelayAAFA
            
            % ======= Get simultaneous time frames =======
            StartIC = round(CleanedMrk(1,2)/(ic_Hdr.samplingfreq/hd_Hdr.samplingfreq));
            StartHD = CleanedMrk(1,3)+DelayAAFA;
            % Time [ms] before marker "1" that can be added to each traces:
            BeforeStart = min(StartIC,StartHD);
            % Delay [ms] between the first and second system starts:
            Delay2syst = max(StartIC,StartHD)-min(StartIC,StartHD);
            
            SimultTime = min(hd_Hdr.nbsamples-StartHD,round(ic_Hdr.nbsamples/(ic_Hdr.samplingfreq/hd_Hdr.samplingfreq))-StartIC);
            
            fprintf('Loading downsampled icEEG data...\n')
            TargetDir = fileparts(strrep(FPOI{n}{2},'1_Raw','2a_Downsampled'));
            load(fullfile(TargetDir,'icEEGi.mat'),'ic_EEGi')
            
            % ======= Align data =======
            hd_EEGaligned = hd_EEG(:,(StartHD-BeforeStart+1):(StartHD+SimultTime-1));
            ic_EEGialigned = ic_EEGi(:,(StartIC-BeforeStart+1):(StartIC+SimultTime-1));
            
            % ======= Align markers as well =======
            ic_T1_aligned = round(ic_T1/(ic_Hdr.samplingfreq/hd_Hdr.samplingfreq))-StartIC+BeforeStart+1;
            ic_T2_aligned = round(ic_T2/(ic_Hdr.samplingfreq/hd_Hdr.samplingfreq))-StartIC+BeforeStart+1;
            
            hd_T1_aligned = hd_T1+DelayAAFA-StartHD+BeforeStart+1;
            hd_T2_aligned = hd_T2+DelayAAFA-StartHD+BeforeStart+1;
            
            % ======= Rename markers to distinguish them from manual and Persyst ones =======
            ic_L_aligned = ic_L;
            for ml = 1:length(ic_L_aligned)
                if ~isempty(regexp(ic_L{ml},'"', 'once'))
                    ic_L_aligned{ml} = ['"EprimeMicromed_',regexprep(ic_L{ml},'"',''),'"'];
                else
                    ic_L_aligned{ml} = ['EprimeMicromed_',ic_L{ml}];
                end
            end
            hd_L_aligned = hd_L;
            for ml = 1:length(hd_L_aligned)
                if ~isempty(regexp(hd_L{ml},'"', 'once'))
                    hd_L_aligned{ml} = ['"EprimeNetStation_',regexprep(hd_L{ml},'"',''),'"'];
                else
                    hd_L_aligned{ml} = ['EprimeNetStation_',hd_L{ml}];
                end
            end
            
            % ===================================================================
            % ======= Keep manual markers added with OLD REALIGN strategy =======
            % ===================================================================
            
            fprintf('Aligning manual markers roughly...\n')
            if exist([strrep(fileparts(FPOI{n}{2}),'1_Raw','2_Resampled'),filesep,'icEEG',filesep,'icEEG.sef.mrk'],'file')==2 && exist([strrep(fileparts(FPOI{n}{1}),'1_Raw','2_Resampled'),filesep,'hdEEG',filesep,'hdEEG.sef.mrk'],'file')==2
                % #RM@FBMlab: does not work anymore because we have moved
                % folder "2_Resampled" from E:\  to I:\, however there are
                % no "manual" (in the sense "clinical") markers for these
                % subjects yet, so it does not matter !
                
                [man_ic_T1,man_ic_T2,man_ic_L] = read_mrk_file([strrep(fileparts(FPOI{n}{2}),'1_Raw','2_Resampled'),filesep,'icEEG',filesep,'icEEG.sef.mrk']);
                [man_hd_T1,man_hd_T2,man_hd_L] = read_mrk_file([strrep(fileparts(FPOI{n}{1}),'1_Raw','2_Resampled'),filesep,'hdEEG',filesep,'hdEEG.sef.mrk']);
                
                % % Check that duration is 0 !
                if ~all((man_ic_T1-man_ic_T2)==0)
                    cprintf('green_','Some manual icEEG markers have duration > 0 !\n')
                    Temp = unique(man_ic_L((man_ic_T1-man_ic_T2)~=0));
                    for asdf = 1:length(Temp)
                        cprintf('cyan',[Temp{asdf},'\n'])
                    end
                end
                if ~all((man_hd_T1-man_hd_T2)==0)
                    cprintf('green_','Some manual hdEEG markers have duration > 0 !\n')
                    Temp = unique(man_hd_L((man_hd_T1-man_hd_T2)~=0));
                    for asdf = 1:length(Temp)
                        cprintf('cyan',[Temp{asdf},'\n'])
                    end
                end
                
                % ======= icEEG =======
                man_ic_T1_aligned = man_ic_T1+BeforeStart;
                man_ic_T2_aligned = man_ic_T2+BeforeStart;
                
                % ======= hdEEG =======
                man_hd_T1_aligned = man_hd_T1+BeforeStart;
                man_hd_T2_aligned = man_hd_T2+BeforeStart;
                
                % ======= Rename markers to distinguish them from E-Prime / MicroMed / NetStation and Persyst ones =======
                
                man_ic_L_aligned = man_ic_L;
                for ml = 1:length(man_ic_L_aligned)
                    if ~isempty(regexp(man_ic_L{ml},'spike cluster*','once'))
                        if ~isempty(regexp(man_ic_L{ml},'"', 'once'))
                            man_ic_L_aligned{ml} = ['"Epitome_',regexprep(man_ic_L{ml},'"',''),'"'];
                        else
                            man_ic_L_aligned{ml} = ['Epitome_',man_ic_L{ml}];
                        end
                    else
                        if ~isempty(regexp(man_ic_L{ml},'"', 'once'))
                            man_ic_L_aligned{ml} = ['"ManSEEG_',regexprep(man_ic_L{ml},'"',''),'"'];
                        else
                            man_ic_L_aligned{ml} = ['ManSEEG_',man_ic_L{ml}];
                        end
                    end
                end
                man_hd_L_aligned = man_hd_L;
                for ml = 1:length(man_hd_L_aligned)
                    if ~isempty(regexp(man_hd_L{ml},'spike cluster*','once'))
                        if ~isempty(regexp(man_hd_L{ml},'"', 'once'))
                            man_hd_L_aligned{ml} = ['"Epitome_',regexprep(man_hd_L{ml},'"',''),'"'];
                        else
                            man_hd_L_aligned{ml} = ['Epitome_',man_hd_L{ml}];
                        end
                    else
                        if ~isempty(regexp(man_hd_L{ml},'"', 'once'))
                            man_hd_L_aligned{ml} = ['"ManScalp_',regexprep(man_hd_L{ml},'"',''),'"'];
                        else
                            man_hd_L_aligned{ml} = ['ManScalp_',man_hd_L{ml}];
                        end
                    end
                end
                
                % =================================
                % ======= MERGE ALL MARKERS =======
                % =================================
                
                % ======= icEEG =======
                ic_T1_aligned_ALL = [ic_T1_aligned;man_ic_T1_aligned];
                ic_T2_aligned_ALL = [ic_T2_aligned;man_ic_T2_aligned];
                ic_L_aligned_ALL = [ic_L_aligned;man_ic_L_aligned];
                
                [ic_T1_aligned_ALL,IdxSorted] = sort(ic_T1_aligned_ALL);
                
                ic_L_aligned_ALL = ic_L_aligned_ALL(IdxSorted);
                ic_T2_aligned_ALL = ic_T2_aligned_ALL(IdxSorted);
                
                % ======= hdEEG =======
                hd_T1_aligned_ALL = [hd_T1_aligned;man_hd_T1_aligned];
                hd_T2_aligned_ALL = [hd_T2_aligned;man_hd_T2_aligned];
                hd_L_aligned_ALL = [hd_L_aligned;man_hd_L_aligned];
                
                [hd_T1_aligned_ALL,IdxSorted] = sort(hd_T1_aligned_ALL);
                
                hd_L_aligned_ALL = hd_L_aligned_ALL(IdxSorted);
                hd_T2_aligned_ALL = hd_T2_aligned_ALL(IdxSorted);
                
            elseif exist([strrep(fileparts(FPOI{n}{2}),'1_Raw','2_Resampled'),filesep,'icEEG',filesep,'icEEG.sef.mrk'],'file')==0 && exist([strrep(fileparts(FPOI{n}{1}),'1_Raw','2_Resampled'),filesep,'hdEEG',filesep,'hdEEG.sef.mrk'],'file')==0
                
                % ======= icEEG =======
                ic_T1_aligned_ALL = ic_T1_aligned;
                ic_T2_aligned_ALL = ic_T2_aligned;
                ic_L_aligned_ALL = ic_L_aligned;
                
                % ======= hdEEG =======
                hd_T1_aligned_ALL = hd_T1_aligned;
                hd_T2_aligned_ALL = hd_T2_aligned;
                hd_L_aligned_ALL = hd_L_aligned;
                
            else
                
                error('Marker files only present in some modalities!')
                
            end
            
            % ============================
            % ======= write output =======
            % ============================
            
            % ======= Matlab format =======
            TargetDir = fileparts(strrep(FPOI{n}{2},'1_Raw','2b_Realigned'));
            mkdir(TargetDir)
            
            [TempFol,Part] = fileparts(TargetDir);
            [TempFol,Type] = fileparts(TempFol);
            [~,Subj] = fileparts(TempFol);
            
            fprintf('Saving icEEG data...\n')
            NewSamplingFreq = hd_EEG_sampling_frequency;
            save(fullfile(TargetDir,['icEEGi_',Subj,'_',Type,'_',Part,'.mat']),'ic_EEGialigned','ic_Hdr','NewSamplingFreq','-v7.3');
            save(fullfile(TargetDir,['icEEGi_',Subj,'_',Type,'_',Part,'_mrk.mat']),'ic_L_aligned_ALL','ic_T1_aligned_ALL','ic_T2_aligned_ALL','-v7.3');
            
            fprintf('Saving hdEEG data...\n')
            save(fullfile(TargetDir,['hdEEG_',Subj,'_',Type,'_',Part,'.mat']),'hd_EEGaligned','hd_Hdr','-v7.3');
            save(fullfile(TargetDir,['hdEEG_',Subj,'_',Type,'_',Part,'_mrk.mat']),'hd_L_aligned_ALL','hd_T1_aligned_ALL','hd_T2_aligned_ALL','-v7.3');
            
            % ======= For Cartool =======
            fprintf('Writing icEEG data for Cartool...\n')
            write_sef(fullfile(TargetDir,['icEEGi_',Subj,'_',Type,'_',Part,'.sef']),ic_EEGialigned',NewSamplingFreq,ic_Hdr.channelnames');
            write_mrk_file_Cartool(fullfile(TargetDir,['icEEGi_',Subj,'_',Type,'_',Part,'.sef.mrk']),ic_T1_aligned_ALL,ic_T2_aligned_ALL,ic_L_aligned_ALL);
            
            fprintf('Writing hdEEG data for Cartool...\n')
            write_sef(fullfile(TargetDir,['hdEEG_',Subj,'_',Type,'_',Part,'.sef']),hd_EEGaligned',hd_Hdr.samplingfreq,hd_Hdr.channelnames');
            write_mrk_file_Cartool(fullfile(TargetDir,['hdEEG_',Subj,'_',Type,'_',Part,'.sef.mrk']),hd_T1_aligned_ALL,hd_T2_aligned_ALL,hd_L_aligned_ALL);
            
            % ======= For AnyWave =======
            
            load('E:\code\MATLAB\MATLAB\dualEEG\hdchs.mat')
            
            fprintf('Writing dual EEG data for AnyWave...\n')
            mat2ades([ic_EEGialigned;hd_EEGaligned],fullfile(TargetDir,['dual_',Subj,'_',Type,'_',Part,]),...
                hd_Hdr.samplingfreq,[ic_Hdr.channelnames';hdchs],...
                [cellstr(repmat('SEEG',length(ic_Hdr.channelnames),1));...
                cellstr(repmat('EEG',length(hdchs),1))]);
            
            Ms = [ic_T1_aligned_ALL;hd_T1_aligned_ALL];
            Me = [ic_T2_aligned_ALL;hd_T2_aligned_ALL];
            ML = [ic_L_aligned_ALL;hd_L_aligned_ALL];
            [~,SortIdx] = sort(Ms);
            write_mrk_file_AnyWave(fullfile(TargetDir,['dual_',Subj,'_',Type,'_',Part,'.ades.mrk'])',...
                ML(SortIdx),Ms(SortIdx)/hd_Hdr.samplingfreq,(Me(SortIdx)-Ms(SortIdx))/hd_Hdr.samplingfreq);
            % marker timings should be in seconds!
            
        else
            TargetDir = fileparts(strrep(FPOI{n}{2},'1_Raw','2b_Realigned'));
            mkdir(TargetDir)
            
            fileID = fopen(fullfile(TargetDir,'error.log'),'w');
            fprintf(fileID,'%s\n','There are no overlapping markers for this recording');
            fclose(fileID);
        end
    end
    
end

%% Check that there are marker files for each modalities (or nothing at all)
for s = 1:length(SubjReadyForCleanCheck)
    
    fprintf('Doing subject %s/%s...\n',num2str(s),num2str(length(SubjReadyForCleanCheck)))
    clear FPOI MergedMrk Header HDeName ICeName hdEEGdur icEEGdur MrkDur Prots
    
    [f,d,fp] = my_recursive_listfiles(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),'');
    Prots = unique(d);
    Prots = Prots(cellfun(@isempty,regexp(Prots,'electrodes_coordinates'))); % ignore "electrodes_coordinates" folder
    Prots = Prots(~strcmp(fullfile(RootDir,SubjDir{SubjReadyForCleanCheck(s)}),Prots)); % if there are files at the subject's root folder level, ignore paths corresponding to these...
    
    Prots = Prots(cellfun(@isempty,regexp(Prots,'DO_NOT_USE'))); % IGNORE ALSO FOLDERS MATCHING THIS PATTERN
    
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
    
    for n = 1:length(Prots)
        
        icMrkFileExist{s}(n) = exist([strrep(fileparts(FPOI{n}{2}),'1_Raw','2_Resampled'),filesep,'icEEG',filesep,'icEEG.sef.mrk'],'file')==2;
        hdMrkFileExist{s}(n) = exist([strrep(fileparts(FPOI{n}{1}),'1_Raw','2_Resampled'),filesep,'hdEEG',filesep,'hdEEG.sef.mrk'],'file')==2;
        
    end
end

for s = 1:length(SubjReadyForCleanCheck)
    
    if isempty(setdiff(1,unique(icMrkFileExist{s}+hdMrkFileExist{s})))
        warning('Problem with subject %s!\n',num2str(s))
    end
    
end


%% Backup and clean marker files (get rid of synchronization triggers, but keep others)

FileList = {'paths_to_mrk_files'};

% Fixed improper read of marker timings (relevant information in newest
% patients is not after "DIN*:" but in "DIN*"):

FileList = {'paths_to_mrk_files_for_patient82_and_patient96'};

for f = 1:length(FileList)
    
    % BKP:
    copyfile(FileList{f},spm_file(FileList{f},'suffix','_BKP_before_marking'));
    
    % get rid of triggers, keep other things:
    [TimeStart, TimeEnd, Label] = read_mrk_Cartool(FileList{f});
    [SyncTriggers,IsNumeric] = dual_read_mkr_NEW(FileList{f});
    
    if ~isempty(TimeStart(~IsNumeric)) % if there remains anything
        write_mrk_file_Cartool(FileList{f},TimeStart(~IsNumeric),TimeEnd(~IsNumeric),Label(~IsNumeric));
    else
        delete(FileList{f});
    end
    
end
