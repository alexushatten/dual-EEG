

cowsay('Fixing ECG signal for patient3...')
subj = 2
% for subj = SubjInd
    
    %% Paths to data
    [~,~,ScalpSefFilepaths] = my_recursive_listfiles(fullfile(RootPath,SubPath{subj},DataType),'hdEEG');
    ScalpSefFilepaths = cellstr(ScalpSefFilepaths);
    ScalpSefFilepaths = ScalpSefFilepaths(~cellfun(@isempty,regexp(ScalpSefFilepaths,'.sef$')));
    SEEGsefFilepaths = regexprep(ScalpSefFilepaths,'hdEEG_','icEEGi_');
    [~,~,XYZfidPath] = my_recursive_listfiles(fullfile(PreprocPath,SubPath{subj},SubPathXYZnorm),'.ToFiducial.xyz');
    
    if length(ScalpSefFilepaths)~=length(SEEGsefFilepaths)
        error('Number of icEEG and hdEEG files differ !')
    end
    
    %% Get data, bandpass & notch filter, downsample and concatenate data
    FileJunctions.start = nan(length(ScalpSefFilepaths),1);
    FileJunctions.end = nan(length(ScalpSefFilepaths),1);
    FileJunctionsNewFreq.start = FileJunctions.start;
    FileJunctionsNewFreq.end = FileJunctions.end;
    Nchans = nan(length(ScalpSefFilepaths),1);
    SamplingFreq = nan(length(ScalpSefFilepaths),1);
    icNchans = nan(length(SEEGsefFilepaths),1);
    icChannelNames = cell(length(SEEGsefFilepaths),1);
    clear MstartF MendF MstartFd MendFd MlabF
    Mstart = []; Mend = []; MstartD = []; MendD = []; Mlab = {};
    Count = 0;
    for f = 1:length(ScalpSefFilepaths)
        %% Get basic file information (header)
        fprintf('================================\n================================\nLoading metadata for:\n"%s"\n================================\n================================\n',ScalpSefFilepaths{f})
        Hdr = dual_load_sef_hdr(ScalpSefFilepaths{f});
        Nchans(f) = Hdr.nbchannels;
        SamplingFreq(f) = Hdr.samplingfreq;
        FileJunctions.start(f) = Count+1;
        FileJunctions.end(f) = Count+Hdr.nbsamples;
        % Get file start and end for downsampled data in advance, because it is
        % easier to do it here:
        %-------------------------------------------------------------------------
        % Intuitively, I thought I would have to use floor() rather than
        % ceil(), but it seems that downsampling output always ends up with 1
        % more sample than the one expected, so ... using ceil() instead /!\
        %-------------------------------------------------------------------------
        FileJunctionsNewFreq.start(f) = ceil(Count/(SamplingFreq(f)/NewFreq))+1;
        FileJunctionsNewFreq.end(f) = ceil(Count/(SamplingFreq(f)/NewFreq))+ceil(Hdr.nbsamples/(SamplingFreq(f)/NewFreq));
        
        icHdr = dual_load_sef_hdr(SEEGsefFilepaths{f});
        if icHdr.nbsamples~=Hdr.nbsamples
            error('Intracranial and scalp EEG do not have the same number of samples !')
        end
        if icHdr.samplingfreq~=Hdr.samplingfreq
            error('Intracranial and scalp EEG do not have the same sampling rate !')
        else
            icSamplingFreq = icHdr.samplingfreq;
        end
        icChannelNames{f,:} = icHdr.channelnames';
        icNchans(f) = icHdr.nbchannels;
        
        %% ==== Markers ====
        %==== ... scalp ====
        [tempstart, tempend, templab] = read_mrk_Cartool([ScalpSefFilepaths{f},'.mrk']);
        %===== ...SEEG =====
        [ictempstart, ictempend, ictemplab] = read_mrk_Cartool([SEEGsefFilepaths{f},'.mrk']);
        
        % Rename labels:
        ictemplab = spm_file(ictemplab,'prefix','SEEG_');
        ictemplab = regexprep(ictemplab,'SEEG_"','"SEEG_');
        templab = spm_file(templab,'prefix','scalp_');
        templab = regexprep(templab,'scalp_"','"scalp_');
        
        % Merge markers:
        tempstart = [tempstart;ictempstart]; %#ok<AGROW>
        tempend = [tempend;ictempend]; %#ok<AGROW>
        templab = [templab;ictemplab]; %#ok<AGROW>
        
        % Sort them:
        [~,IdxSort] = sort(tempstart);
        tempstart = tempstart(IdxSort);
        tempend = tempend(IdxSort);
        templab = templab(IdxSort);
        
        % Adjust marker timings: remember that the rule is "mrk + 1" !
        % (because indexing starts at 1 in Matlab vs 0 in Cartool!
        MstartF{f} = tempstart+1; %#ok<SAGROW>
        MendF{f} = tempend+1; %#ok<SAGROW>
        MlabF{f} = templab; %#ok<SAGROW>
        if ~all(length(tempstart) == [length(tempend),length(templab)])
            % This should not be the case by construction (see
            % read_mrk_Cartool.m), but just to be sure:
            error('Error reading markers: start, end and label vectors do not match !')
        end
        IndicesForThisIter = (length(Mstart)+1):(length(Mstart)+length(tempstart));
        Mstart(IndicesForThisIter,1) = tempstart+Count+1; %#ok<SAGROW>
        Mend(IndicesForThisIter,1) = tempend+Count+1; %#ok<SAGROW>
        Mlab(IndicesForThisIter,1) = templab; %#ok<SAGROW>
        
        %% ==== Downsample marker timings ====
        % This should be done in advance, because it is easier to do it here:
        % ... Easy part, because it restarts each time (although these variables
        % should in principle not used...):
        MstartFd{f} = ceil(MstartF{f}-1/(SamplingFreq(f)/NewFreq))+1; %#ok<SAGROW>
        MendFd{f} = ceil(MendF{f}-1/(SamplingFreq(f)/NewFreq))+1; %#ok<SAGROW>
        % ... More difficult one, because we might have lost some samples in the
        % process:
        MstartD(IndicesForThisIter,1) = ceil(tempstart/(SamplingFreq(f)/NewFreq))+ceil(Count/(SamplingFreq(f)/NewFreq))+1; %#ok<SAGROW>
        MendD(IndicesForThisIter,1) = ceil(tempend/(SamplingFreq(f)/NewFreq))+ceil(Count/(SamplingFreq(f)/NewFreq))+1; %#ok<SAGROW>
        
        % update Count
        Count = Count+Hdr.nbsamples;
    end
    %% Check for concatenation issues
    FileDurations = FileJunctions.end-FileJunctions.start+1;
    FileDurationsNewFreq = FileJunctionsNewFreq.end-FileJunctionsNewFreq.start+1;
    if ~all(Nchans(1)==Nchans)
        error('Files do not all have the same number of channels !')
    else
        Nchans = Nchans(1);
        if ~all(FileJunctions.start(2:end)-FileJunctions.end(1:end-1)) || ~all(FileJunctionsNewFreq.start(2:end)-FileJunctionsNewFreq.end(1:end-1))
            error('File starts and ends do not match, check for errors !')
        end
        if ~all(SamplingFreq(1)==SamplingFreq)
            error('Sampling rate is not the same across files !')
        else
            SamplingFreq = SamplingFreq(1);
            if (rem((SamplingFreq/NewFreq),1)>0)
                error('Targeted downsampling frequency is not an integer!')
            else
                % the downsampling is done on each file separately, so so we might loose more samples than expected !
                EEGd = nan(Nchans,FileJunctionsNewFreq.end(end));
                ECGd = nan(1,FileJunctionsNewFreq.end(end));
                ECGdAbsFilt = nan(1,FileJunctionsNewFreq.end(end));
                if ~all(icNchans(1)==icNchans)
                    for f = 1:length(SEEGsefFilepaths)
                        %                     % filter out "none" channel
                        %                     icChannelNames{f} = icChannelNames{f}(cellfun(@isempty,regexp(icChannelNames{f},'^none$')));
                        %                     % filter out "ecg" channel, because it might be not
                        %                     % present in some files:
                        %                     icChannelNames{f} = icChannelNames{f}(cellfun(@isempty,regexpi(icChannelNames{f},'^ecg$')));
                        [Bipoles,ElecMatch,labelsB,labelsB2] = bipolar_montage(icChannelNames{f},2);
                        SuspiciousChannels = icChannelNames{f}(cellfun(@isempty,match_vectors(icChannelNames{f},unique(labelsB2(:)),1)));
                        SuspiciousChannelsIdx = match_vectors(SuspiciousChannels,icChannelNames{f},1);
                        icChannelNames{f} = icChannelNames{f}(setdiff(1:length(icChannelNames{f}),SuspiciousChannelsIdx));
                    end
                    icNchans = cellfun(@length,icChannelNames);
                    if ~all(icNchans(1)==icNchans)
                        error('Number of channels vary across intracranial files, cannot pre-allocate !')
                    end
                end
                icNchans = icNchans(1);
                icChannelNames = icChannelNames{1};
                icEEGd = nan(icNchans,FileJunctionsNewFreq.end(end));
                %             EEGd = nan(Nchans,sum(FileDurationsNewFreq));
                %             EEGd = nan(Nchans,sum(floor(FileDurations/(SamplingFreq/NewFreq))));
            end
            EEG = nan(Nchans,FileJunctions.end(end));
            ECG = nan(1,FileJunctions.end(end));
            ECGabsFilt = nan(1,FileJunctions.end(end));
            icEEG = nan(icNchans,FileJunctions.end(end));
        end
    end
    % Extract cardiac data if available, filter, downsample and concatenate:
    CardiacChanSplit = false;
    for f = 1:length(ScalpSefFilepaths)
        fprintf('================================\n================================\nLoading data part %d / %d...\n================================\n================================\n',f,length(ScalpSefFilepaths))
        
        fprintf('Loading scalp data in "%s"...\n',ScalpSefFilepaths{f})
        % load raw hdEEG traces
        [~,Hdr] = dual_load_sef(ScalpSefFilepaths{f});
        
        % Get also channel names
        ChannelNames = Hdr.channelnames;
        
        % load also raw icEEG traces (icEEG and hdEEG traces
        % SHOULD ALREADY BE ALIGNED !):
        fprintf('Loading SEEG data in "%s"...\n',SEEGsefFilepaths{f})
        [icData,icHdr] = dual_load_sef(SEEGsefFilepaths{f});
        
        icChannelNamesAll = icHdr.channelnames';
        
        %% Find cardiac channel in icEEG if any
        % ... will be useful for ICA later:
        
        % check if there is a channel that cannot be paired with others
        [Bipoles,ElecMatch,labelsB,labelsB2] = bipolar_montage(icChannelNamesAll,2);
        try
            SuspiciousChannels = icChannelNamesAll(cellfun(@isempty,match_vectors(icChannelNamesAll,unique(labelsB2(:)),1)));
            SuspiciousChannelsIdx = match_vectors(SuspiciousChannels,icChannelNamesAll,1);
%             CardiacIdx = ~cellfun(@isempty,regexpi(SuspiciousChannels,'^ecg( )*$'));
            CardiacIdx = ~cellfun(@isempty,regexpi(SuspiciousChannels,'^ecg( |\d)*$')); % in one subject (RD), the ECG channel is labeled "ECG1"
            CardiacIdx = SuspiciousChannelsIdx(CardiacIdx);
        catch
            SuspiciousChannels = [];
            SuspiciousChannelsIdx = [];
            CardiacIdx = [];
        end
        
        % Sometimes the ECG channel has no name and is just labeled "e" +
        % electrode number:
        if isempty(CardiacIdx) && ~isempty(SuspiciousChannels)
            CardiacIdx = ~cellfun(@isempty,regexpi(SuspiciousChannels,'^e\d+'));
            warning('Channel "%s" identified as ECG channel, please make sure it is correct !',SuspiciousChannels{CardiacIdx})
            %         % Check if marker channel (should not be the case in principle):
            %         regexpi('^mkr')
            CardiacIdx = SuspiciousChannelsIdx(CardiacIdx);
        end
        
        % Sometimes the ECG channel is split in positive / negative parts:
        if isempty(CardiacIdx) && ~isempty(SuspiciousChannels)
            CardiacIdx = ~cellfun(@isempty,regexpi(SuspiciousChannels,'^ecg( )*(+|-)( )*$'));
            warning('Channel "%s" identified as ECG channel, please make sure it is correct !',SuspiciousChannels{CardiacIdx})
            %         % Check if marker channel (should not be the case in principle):
            %         regexpi('^mkr')
            CardiacIdx = SuspiciousChannelsIdx(CardiacIdx);
            if ~isempty(CardiacIdx)
                if length(CardiacIdx)==2
                    CardiacChanSplit = true;
                    warning('ECG channel found but split in +/-, will calculate difference between two channels...')
                elseif length(CardiacIdx)>2
                    error('Multiple channels match regexp, check file %s and code !',SEEGsefFilepaths{f})
                end
            else
                warning('Could not find ECG channel for file %s !',SEEGsefFilepaths{f});
            end
        end
        
        %% Filtering (bandpass + notch) hdEEG, icEEG & cardiac
        fprintf('================================\n================================\nFiltering part %d / %d...\n================================\n================================\n',f,length(ScalpSefFilepaths))
        fprintf('Butterworth bandpass filter order 4, highpass = %d, lowpass = %d\n',HighPassBand,LowPassBand)
%         EEG(:,FileJunctions.start(f):FileJunctions.end(f)) = ...
%             dual_filt(Data,SamplingFreq,[HighPassBand LowPassBand],LineNoise);
        % if any ECG, filter it as well
        if ~isempty(CardiacIdx)
            if CardiacChanSplit % if ECG is split in +/-, make difference (bipolar montage) to get ECG
                ECG(FileJunctions.start(f):FileJunctions.end(f)) = ...
                    dual_filt(icData(CardiacIdx(1),:)-icData(CardiacIdx(2),:),SamplingFreq,[HighPassBand LowPassBand],LineNoise);
            else
                ECG(FileJunctions.start(f):FileJunctions.end(f)) = ...
                    dual_filt(icData(CardiacIdx,:),SamplingFreq,[HighPassBand LowPassBand],LineNoise);
            end
            ECGabsFilt(FileJunctions.start(f):FileJunctions.end(f)) = ...
                ft_preproc_bandpassfilter(abs(ECG(FileJunctions.start(f):FileJunctions.end(f))),SamplingFreq,[ECGhp ECGlp],ECGbpOrder,'but','twopass');
        end
        % filter SEEG as well
%         icEEG(:,FileJunctions.start(f):FileJunctions.end(f)) = ...
%             dual_filt(icData(1:length(icChannelNames),:),SamplingFreq,[HighPassBand LowPassBand],LineNoise);
        
        %% Downsampling filtered hdEEG, icEEG & cardiac
        % Done here, because:
        % - speeds up ICA and still allows to reconstruct at higher sampling rate
        % based on ICA weights (should be fine as long as sampling rate after downsampling is
        % above 120 Hz (see https://neuroimage.usc.edu/forums/t/downsampling-data-before-ica-processing/3064)
        % - should be done before concatenation because of preliminary low-pass filtering
        
        fprintf('================================\n================================\nDownsampling part %d / %d to %d Hz\nusing Cascaded-Integrator Comb filters...\n================================\n================================\n',f,length(ScalpSefFilepaths),NewFreq)
%         EEGd(:,FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
%             downsample_CIC_EEG( EEG(:,FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
        if ~isempty(CardiacIdx)
            ECGd(FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
                downsample_CIC_EEG( ECG(FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
            ECGdAbsFilt(FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
                downsample_CIC_EEG( ECGabsFilt(FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
        end
%         icEEGd(:,FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
%             downsample_CIC_EEG( icEEG(:,FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
        
    end
    
    %% Export in Matlab & Cartool format, although the latter might not be able to open all of them (2GB limit)...
    fprintf('================================\n================================\nSaving output...\n================================\n================================\n')
    
    OutFilePath = regexprep(spm_file(spm_file(ScalpSefFilepaths{1},'path'),'path'),'2b_Realigned','3_Preprocessed');
%     mkdir(OutFilePath)
    %============ save scalp ... =============
%     fprintf('Saving scalp...\n')
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    % filtered EEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise)];
    load([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEG','SamplingFreq','ChannelNames','FileJunctions','Nchans',...
        'MstartF','MendF','MlabF','Mstart','Mend','Mlab');
    % filtered and downsampled EEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    load([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGd','NewFreq','ChannelNames','FileJunctionsNewFreq','Nchans',...
        'MstartFd','MendFd','MlabF','MstartD','MendD','Mlab');
    
    %============= ... and separate SEEG =============
%     fprintf('Saving SEEG...\n')
    OutFilename = spm_file(SEEGsefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    % filtered SEEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise)];
    load([OutFilePath,filesep,OutFilename,'.mat'],...
        'icEEG','SamplingFreq','icChannelNames','FileJunctions','Nchans',...
        'MstartF','MendF','MlabF','Mstart','Mend','Mlab');

    % filtered and downsampled SEEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    load([OutFilePath,filesep,OutFilename,'.mat'],...
        'icEEGd','NewFreq','icChannelNames','FileJunctionsNewFreq','Nchans',...
        'MstartFd','MendFd','MlabF','MstartD','MendD','Mlab');

    %===== All waves for AnyWave =====
    EGI = get_EGI_257ch_sensors_info;
    % load('E:\code\MATLAB\MATLAB\dualEEG\hdchs.mat') % replaced with get_EGI_257ch_sensors_info.m and look for field
    
    fprintf('Writing dual file for AnyWave...\n')
    % filtered
    OutFilename = ['dual_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise)];
    if ~all(isempty(ECG))
        
        ECGnanZeroed = ECG;
        ECGnanZeroed(isnan(ECGnanZeroed))=0;
        
        mat2ades([icEEG;ECGnanZeroed;EEG],fullfile(OutFilePath,OutFilename),...
            NewFreq,[icChannelNames;{'ECG'};EGI.clinicalname],...
            [cellstr(repmat('SEEG',length(icChannelNames),1));...
            {'ECG'};...
            cellstr(repmat('EEG',length(EGI.clinicalname),1))]);
    else
        mat2ades([icEEG;EEG],fullfile(OutFilePath,OutFilename),...
            NewFreq,[icChannelNames;EGI.clinicalname],...
            [cellstr(repmat('SEEG',length(icChannelNames),1));...
            cellstr(repmat('EEG',length(EGI.clinicalname),1))]);
    end
    write_mrk_file_AnyWave([fullfile(OutFilePath,OutFilename),'.ades.mrk'],...
        [Mlab;repmat({'"File_junction"'},length(FileJunctions.start)-1,1)],...
        [Mstart-1;FileJunctions.start(2:end)-1],...
        [Mend-1;FileJunctions.start(2:end)-1]-[Mstart-1;FileJunctions.start(2:end)-1]);
    
    % filtered and downsampled
    OutFilename = ['dual_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    if ~all(isempty(ECG))
        
        ECGdNanZeroed = ECGd;
        ECGdNanZeroed(isnan(ECGdNanZeroed))=0;
        
        mat2ades([icEEGd;ECGdNanZeroed;EEGd],fullfile(OutFilePath,OutFilename),...
            NewFreq,[icChannelNames;{'ECG'};EGI.clinicalname],...
            [cellstr(repmat('SEEG',length(icChannelNames),1));...
            {'ECG'};...
            cellstr(repmat('EEG',length(EGI.clinicalname),1))]);
    else
        mat2ades([icEEGd;EEGd],fullfile(OutFilePath,OutFilename),...
            NewFreq,[icChannelNames;EGI.clinicalname],...
            [cellstr(repmat('SEEG',length(icChannelNames),1));...
            cellstr(repmat('EEG',length(EGI.clinicalname),1))]);
    end
    write_mrk_file_AnyWave([fullfile(OutFilePath,OutFilename),'.ades.mrk'],...
        [Mlab;repmat({'"File_junction"'},length(FileJunctionsNewFreq.start)-1,1)],...
        [MstartD-1;FileJunctionsNewFreq.start(2:end)-1],...
        [MendD-1;FileJunctionsNewFreq.start(2:end)-1]-[MstartD-1;FileJunctionsNewFreq.start(2:end)-1]);
    
    % free some memory:
    clear EEG
    
    % ===== save also separate cardiac signal... =====
    if ~all(isempty(ECG))
        fprintf('Saving ECG...\n')
        % filtered ECG:
        OutFilename = ['cardiac_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise)];
        save([OutFilePath,filesep,OutFilename,'.mat'],...
            'ECG','SamplingFreq','FileJunctions','-v7.3');
        ECGnanZeroed = ECG;
        ECGnanZeroed(isnan(ECGnanZeroed))=0;
        write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
            ECGnanZeroed',SamplingFreq,{'ECG'});
        % filtered and downsampled ECG:
        OutFilename = ['cardiac_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
        save([OutFilePath,filesep,OutFilename,'.mat'],...
            'ECGd','NewFreq','FileJunctionsNewFreq','-v7.3');
        ECGdNanZeroed = ECGd;
        ECGdNanZeroed(isnan(ECGdNanZeroed))=0;
        write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
            ECGdNanZeroed',NewFreq,{'ECG'});
        
        % filtered ECG:
        OutFilename = ['cardiac_',TempOutFilename{2},'_',TempOutFilename{3},'_bp1-15hz'];
        save([OutFilePath,filesep,OutFilename,'.mat'],...
            'ECGabsFilt','SamplingFreq','FileJunctions','-v7.3');
        ECGnanZeroed = ECGabsFilt;
        ECGnanZeroed(isnan(ECGnanZeroed))=0;
        write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
            ECGnanZeroed',SamplingFreq,{'ECG'});
        % filtered and downsampled ECG:
        OutFilename = ['cardiac_',TempOutFilename{2},'_',TempOutFilename{3},'_bp1-15hz_down',num2str(NewFreq),'hz'];
        save([OutFilePath,filesep,OutFilename,'.mat'],...
            'ECGdAbsFilt','NewFreq','FileJunctionsNewFreq','-v7.3');
        ECGdNanZeroed = ECGdAbsFilt;
        ECGdNanZeroed(isnan(ECGdNanZeroed))=0;
        write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
            ECGdNanZeroed',NewFreq,{'ECG'});
    end
    
    % free some memory:
    clear ECG

    % free some memory:
    clear icEEG
    
    % Save results
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_artefacts'];
    
    fprintf('================================\n================================\nArtefacts detection done, saving...\n================================\n================================\n')
    load(fullfile(OutFilePath,[OutFilename,'.mat']),...
        'Bad_Channels', 'Good_Channels', 'Bad_Segments', 'Good_Segments',...
        'EGI', 'Metrics', 'ArtDetectResults', 'OnsetsBadIdx', 'OffsetsBadIdx')
    
    %% ICA
    
%     fprintf('================================\n================================\nStarting ICA decomposition...\n================================\n================================\n\n')
    %=========== ICA decomposition ===========
%     [W,Winv,Activations,Topos,CompWithCardiacFreq,FreqPower,FreqList] = ...
%         my_ICA( EEGd(Good_Channels,Good_Segments), NewFreq );
    
%     fprintf('================================\n================================\nICA decomposition done, saving...\n================================\n================================\n')
    % Save results
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_ICA'];
    
    load(fullfile(OutFilePath,[OutFilename,'.mat']),'W','Winv','Activations','Topos','CompWithCardiacFreq','FreqPower','FreqList')
    
    % Adjust markers after bad segments excision:
    [AdjMstartD,AdjMendD,AdjMlab] = adjust_mrk_4_bad_segments(MstartD,MendD,Mlab,OnsetsBadIdx,OffsetsBadIdx);
    
    % Write XYZ without bad channels for visualization in Cartool:
%     write_xyz_file([OutFilePath,filesep,'Template_cap_without_bad_electrodes.xyz'],EGI.X(Good_Channels),EGI.Y(Good_Channels),EGI.Z(Good_Channels),EGI.clinicalname(Good_Channels))
    
    % Write activations time course and topography for visualization in Cartool:
    CompNames = regexprep(cellstr([repmat('C',size(Activations,1),1),num2str((1:size(Activations,1))')]),' ','');
%     write_sef(fullfile(OutFilePath,[OutFilename,'_time_courses.sef']),Activations',NewFreq,CompNames);
    % Below we add zeros for the first "time frame" because Cartool starts
    % at 0. Otherwise component 1 is at time frame 0, component 2 at time
    % frame 1, ... and so on... This avoids making stupid errors when
    % looking at IC topography:
%     write_sef(fullfile(OutFilePath,[OutFilename,'_topo.sef']),[zeros(1,size(Topos,2));Topos],NewFreq,ChannelNames(:));
%     write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'_time_courses.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % to directly open ICs topography with appropriate XYZ:
%     write_lm_file(fullfile(OutFilePath,[OutFilename,'_topo.lm']),[{fullfile(OutFilePath,[OutFilename,'_topo.sef'])};{[OutFilePath,filesep,'Template_cap_without_bad_electrodes.xyz']}]);
    
    %========= Check correlation with bipolar montage of frontal electrode =========
    %====================== for ocular and muscular artefacts ======================
    % Construct EOG and EMG channels
    Chan4emg = {'FP2'	'F8'
        'F8'	'T8'
        'FP1'	'F7'
        'F7'	'T7'};
    Lchan4emg = match_vectors(Chan4emg(:,1),EGI.clinicalname,1);
    Rchan4emg = match_vectors(Chan4emg(:,2),EGI.clinicalname,1);
    if ~isempty(intersect(Lchan4emg,Bad_Channels)) || ~isempty(intersect(Rchan4emg,Bad_Channels))
        warning('Some of the frontotemporal EEG channels used for bipolar montage highlighting muscular artefacts were marked as bad!')
    end
    EMGlike = EEGd(Lchan4emg,Good_Segments)-EEGd(Rchan4emg,Good_Segments);
    Chan4eog = {'F8'	'T8'
        'FP2'	'F10'
        'FP2'	'F8'
        'FP2'	'F4'
        'F7'	'T7'
        'FP1'	'F3'
        'FP1'	'F7'
        'FP1'	'F9'};
    Lchan4eog = match_vectors(Chan4eog(:,1),EGI.clinicalname,1);
    Rchan4eog = match_vectors(Chan4eog(:,2),EGI.clinicalname,1);
    if ~isempty(intersect(Lchan4eog,Bad_Channels)) || ~isempty(intersect(Rchan4eog,Bad_Channels))
        warning('Some of the frontal EEG channels used for bipolar montage highlighting ocular artefacts were marked as bad!')
    end
    EOGlike = EEGd(Lchan4eog,Good_Segments)-EEGd(Rchan4eog,Good_Segments);
    
    CompEOGcorr = corr(Activations',EOGlike');
    [~,CompEOGcorrMidx] = sort(abs(mean(CompEOGcorr,2)),'descend'); % could also be mean(abs()) but a given IC should in principle be correlated the same way with all bipolar channels because they share very similar info
    CompWithEOG = [EOGlike;Activations(CompEOGcorrMidx,:)]; % Good_Segments not needed here, bad segments already excised
    CompWithEOGlab = [regexprep(cellstr([char(Chan4eog(:,1)),repmat('-',size(Chan4eog,1),1),char(Chan4eog(:,2))]),' ','');CompNames(CompEOGcorrMidx)];
%     write_sef(fullfile(OutFilePath,[OutFilename,'_vs_EOG.sef']),CompWithEOG',NewFreq,CompWithEOGlab);
%     write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'_vs_EOG.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    CompEMGcorr = corr(Activations',EMGlike');
    [~,CompEMGcorrMidx] = sort(abs(mean(CompEMGcorr,2)),'descend'); % could also be mean(abs()) but a given IC should in principle be correlated the same way with all bipolar channels because they share very similar info
    CompWithEMG = [EMGlike;Activations(CompEMGcorrMidx,:)]; % Good_Segments not needed here, bad segments already excised
    CompWithEMGlab = [regexprep(cellstr([char(Chan4emg(:,1)),repmat('-',size(Chan4emg,1),1),char(Chan4emg(:,2))]),' ','');CompNames(CompEMGcorrMidx)];
%     write_sef(fullfile(OutFilePath,[OutFilename,'_vs_EMG.sef']),CompWithEMG',NewFreq,CompWithEMGlab);
%     write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'_vs_EMG.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    ICAcorr.CompEOGcorr = CompEOGcorr;
    ICAcorr.CompEMGcorr = CompEMGcorr;
    
    %========= Check correlation with cardiac signal from SEEG =========
    % ... if available:
    if ~all(isnan(ECGdAbsFilt))
        ECGsig = ECGdAbsFilt(Good_Segments);
        % if NaN then set to 0:
        ECGsig(isnan(ECGsig))=0;
        % Take segments corresponding to ICA activations and look at
        % correlation with ICA time courses:
        CompECGcorr = corr(Activations',ECGsig');
        [~,CompECGcorrMidx] = sort(abs(CompECGcorr));
        CompWithECG = [ECGsig;Activations(CompECGcorrMidx,:)];
        CompWithECGlab = [{'ECG'};CompNames(CompECGcorrMidx)];
        write_sef(fullfile(OutFilePath,[OutFilename,'_vs_ECG.sef']),CompWithECG',NewFreq,CompWithECGlab);
        write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'_vs_ECG.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
        ICAcorr.CompECGcorr = CompECGcorr;
    else
        warning('No ECG channel was found for this subject... Cannot look at correlation between ICA components and cardiac signal.')
    end
    
    %========= FASTER metrics on ICs =========
    [ MeanGradientHF, SpatialKurtosis, IC_hurst_exp, Zscores ] = ICA_FASTER_metrics( Activations, Winv );
    
    %========= For ICs review =========
    ICAinfos = [];
    ICAinfos.FASTER.MeanGradientHF = MeanGradientHF;
    ICAinfos.FASTER.SpatialKurtosis = SpatialKurtosis;
    ICAinfos.FASTER.IC_hurst_exp = IC_hurst_exp;
    ICAinfos.FASTER.Zscores = Zscores;
    ICAinfos.SamplingFreq = NewFreq;
    ICAinfos.ICAcorr = ICAcorr;
    ICAinfos.TopoPath = fullfile(OutFilePath,[OutFilename,'_topo.lm']);
    ICAinfos.ActivationPath = fullfile(OutFilePath,[OutFilename,'_time_courses.sef']);
    ICAinfos.EOGlikePath = fullfile(OutFilePath,[OutFilename,'_vs_EOG.sef']);
    ICAinfos.EMGlikePath = fullfile(OutFilePath,[OutFilename,'_vs_EMG.sef']);
    if ~all(isnan(ECGdAbsFilt))
        ICAinfos.LPFabsECGpath = fullfile(OutFilePath,[OutFilename,'_vs_ECG.sef']);
    end
    
    save(fullfile(OutFilePath,[OutFilename,'_infos.mat']),'ICAinfos');
    
    %========= Write output in AnyWave format with only good segments and good channels =========
    OutFilename = ['dual_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz_good_hdEEG_segments_and_channels'];
    if ~all(isnan(ECGdAbsFilt))
        
        ECGdNanZeroed = ECGd;
        ECGdNanZeroed(isnan(ECGdNanZeroed))=0;
        
        mat2ades([icEEGd(:,Good_Segments);ECGd(:,Good_Segments);EEGd(Good_Channels,Good_Segments)],...
            fullfile(OutFilePath,OutFilename),...
            NewFreq,[icChannelNames;{'ECG'};EGI.clinicalname(Good_Channels)],...
            [cellstr(repmat('SEEG',length(icChannelNames),1));...
            {'ECG'};...
            cellstr(repmat('EEG',length(EGI.clinicalname(Good_Channels)),1))]);
    else
        mat2ades([icEEGd(:,Good_Segments);EEGd(Good_Channels,Good_Segments)],...
            fullfile(OutFilePath,OutFilename),...
            NewFreq,[icChannelNames;EGI.clinicalname(Good_Channels)],...
            [cellstr(repmat('SEEG',length(icChannelNames),1));...
            cellstr(repmat('EEG',length(EGI.clinicalname(Good_Channels)),1))]);
    end
    write_mrk_file_AnyWave([fullfile(OutFilePath,OutFilename),'.ades.mrk'],...
        AdjMlab,AdjMstartD-1,(AdjMendD-1)-(AdjMstartD-1));
    
    
    warning('ECG in patient3 fixed !...')
    
% end