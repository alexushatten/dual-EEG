
% Preprocessing pipeline for dual EEG
clear; clc

%% Parameters
RefChan = 257;
NewFreq = 200;
HighPassBand = 1;
LowPassBand = 45;
LineNoise = 50;
ECGhp = 1;
ECGlp = 15;
ECGbpOrder = 2;
% ... SHOULD THE FILTERS BE DIFFERENT FOR SEEG AND ECG?? EVEN IF A
% HIGH-PASS FILTER IS APPLIED AT 1Hz, THE DEFLECTIONS ELICITED BY HEART
% BEATS ARE STILL CLEARLY VISIBLE... THAT MIGHT BE BECAUSE OF THE FILTER
% ORDER (4), BECAUSE USING A 2ND ORDER BUTTERWORTH WILL DRASTICALLY CHANGE
% RESULTS (THIS IS WHAT CARTOOL DOES !!)...

RootPath = 'E:\dualEEG\2b_Realigned';
DataType = 'resting_state';
Recordings2exclude = 'DO_NOT_USE';
SubPath = {'patient50'
    'patient48'
    'patient82'
    'patient96'};
PreprocPath = 'E:\dualEEG\3_Preprocessed';
SubPathXYZnorm = 'hdEEG_cap_coregistered_fiducials';
SubjInd = 1:length(SubPath);

%% Part 0 (create output directories to store "normalized" (to fiducial) head caps
for subj = SubjInd
    mkdir([PreprocPath,filesep,SubPath{subj},filesep,SubPathXYZnorm])
end
cprintf('cyan','Now interpolate dummy EEG traces in Cartool\nwithout cleaning temporary files and move\n***.FromFiducial.xyz and ***.ToFiducial.xyz files\nto hdEEG_cap_coregistered_fiducials folder!\n')

%% Part I (... until ICA components review and selection)

diary(fullfile(PreprocPath,filesep,'Log_pipeline_v1.txt'))
cowsay('Starting EEG processing pipeline V1...')
for subj = SubjInd
    
    fprintf('%s\n',datestr(now,'yyyy.mm.dd, hh:MM:ss','local'))
    fprintf('Doing subject %d/%d...\n',subj,length(SubjInd))
    %% Paths to data
    [~,~,ScalpSefFilepaths] = my_recursive_listfiles(fullfile(RootPath,SubPath{subj},DataType),'hdEEG');    
    ScalpSefFilepaths = cellstr(ScalpSefFilepaths);
    ScalpSefFilepaths = ScalpSefFilepaths(~cellfun(@isempty,regexp(ScalpSefFilepaths,'.sef$')));
    
    % Filter out recordings that should not be analyzed:
    ScalpSefFilepaths = ScalpSefFilepaths(cellfun(@isempty,regexp(ScalpSefFilepaths,Recordings2exclude)));
    
    SEEGsefFilepaths = regexprep(ScalpSefFilepaths,'hdEEG_','icEEGi_');
    
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
    fprintf('================================\n================================\nLoading metadata infos for files...\n================================\n================================\n')
    for f = 1:length(ScalpSefFilepaths)
        %% Get basic file information (header)
        
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
        [Data,Hdr] = dual_load_sef(ScalpSefFilepaths{f});
        
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
            if sum(CardiacIdx)>0
                warning('Channel "%s" identified as ECG channel, please make sure it is correct !',SuspiciousChannels{CardiacIdx})
            end
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
        EEG(:,FileJunctions.start(f):FileJunctions.end(f)) = ...
            dual_filt(Data,SamplingFreq,[HighPassBand LowPassBand],LineNoise);
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
        icEEG(:,FileJunctions.start(f):FileJunctions.end(f)) = ...
            dual_filt(icData(1:length(icChannelNames),:),SamplingFreq,[HighPassBand LowPassBand],LineNoise);
        
        %% Downsampling filtered hdEEG, icEEG & cardiac
        % Done here, because:
        % - speeds up ICA and still allows to reconstruct at higher sampling rate
        % based on ICA weights (should be fine as long as sampling rate after downsampling is
        % above 120 Hz (see https://neuroimage.usc.edu/forums/t/downsampling-data-before-ica-processing/3064)
        % - should be done before concatenation because of preliminary low-pass filtering
        
        fprintf('================================\n================================\nDownsampling part %d / %d to %d Hz\nusing Cascaded-Integrator Comb filters...\n================================\n================================\n',f,length(ScalpSefFilepaths),NewFreq)
        EEGd(:,FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
            downsample_CIC_EEG( EEG(:,FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
        if ~isempty(CardiacIdx)
            ECGd(FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
                downsample_CIC_EEG( ECG(FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
            ECGdAbsFilt(FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
                downsample_CIC_EEG( ECGabsFilt(FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
        end
        icEEGd(:,FileJunctionsNewFreq.start(f):FileJunctionsNewFreq.end(f)) = ...
            downsample_CIC_EEG( icEEG(:,FileJunctions.start(f):FileJunctions.end(f)), SamplingFreq, NewFreq );
        
    end
    
    %% Export in Matlab & Cartool format, although the latter might not be able to open all of them (2GB limit)...
    fprintf('================================\n================================\nSaving output...\n================================\n================================\n')
    
    OutFilePath = regexprep(spm_file(spm_file(ScalpSefFilepaths{1},'path'),'path'),'2b_Realigned','3_Preprocessed');
    mkdir(OutFilePath)
    %============ save scalp ... =============
    fprintf('Saving scalp...\n')
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    % filtered EEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise)];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEG','SamplingFreq','ChannelNames','FileJunctions','Nchans',...
        'MstartF','MendF','MlabF','Mstart','Mend','Mlab','-v7.3');
    % filtered and downsampled EEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGd','NewFreq','ChannelNames','FileJunctionsNewFreq','Nchans',...
        'MstartFd','MendFd','MlabF','MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGd',NewFreq,ChannelNames);
    % downsampled markers with "file junctions" markers (no need to re-sort them just for Cartool... it will figure it out anyway), just for visualization purposes:
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],...
        [MstartD-1;FileJunctionsNewFreq.start(2:end)-1],[MendD-1;FileJunctionsNewFreq.start(2:end)-1],[Mlab;repmat({'"File_junction"'},length(FileJunctionsNewFreq.start)-1,1)]);
    
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
    
    %============= ... and separate SEEG =============
    fprintf('Saving SEEG...\n')
    OutFilename = spm_file(SEEGsefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    % filtered SEEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise)];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'icEEG','SamplingFreq','icChannelNames','FileJunctions','Nchans',...
        'MstartF','MendF','MlabF','Mstart','Mend','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        icEEG',SamplingFreq,icChannelNames);
    % markers with "file junctions" markers (no need to re-sort them just for Cartool... it will figure it out anyway), just for visualization purposes:
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],...
        [Mstart-1;FileJunctions.start(2:end)-1],[Mend-1;FileJunctions.start(2:end)-1],[Mlab;repmat({'"File_junction"'},length(FileJunctions.start)-1,1)]);
    
    % filtered and downsampled SEEG:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'icEEGd','NewFreq','icChannelNames','FileJunctionsNewFreq','Nchans',...
        'MstartFd','MendFd','MlabF','MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        icEEGd',NewFreq,icChannelNames);
    % downsampled markers with "file junctions" markers (no need to re-sort them just for Cartool... it will figure it out anyway), just for visualization purposes:
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],...
        [MstartD-1;FileJunctionsNewFreq.start(2:end)-1],[MendD-1;FileJunctionsNewFreq.start(2:end)-1],[Mlab;repmat({'"File_junction"'},length(FileJunctionsNewFreq.start)-1,1)]);
    
    % free some memory:
    clear icEEG
    
    %% Detect & remove bad channels & segments using recursive FASTER
    fprintf('================================\n================================\nStarting detection of bad channels and segments...\n================================\n================================\n')
    [Bad_Channels, Good_Channels,...
        Bad_Segments, Good_Segments,...
        Metrics, ArtDetectResults, ArtDetectSettings] = ...
        rm_artefacts(EEGd, NewFreq, RefChan, EGI);
    
    % Display artefact detection settings:
    unfold_struct(ArtDetectSettings)
    
    if ArtDetectResults.Channels.PercentBad>10
        warning('More than 10% of channels are bad, check bad channels / segments detection results!')
    end
    if ArtDetectResults.Segments.PercentBadTF>0.3
        warning('More than 30% of segments are bad, check bad channels / segments detection results!')
    end
    
    %===== Add markers in .mrk (for visual check) =====
    OnsetsBadIdx = find(diff(Bad_Segments)>0);
    OffsetsBadIdx = find(diff(Bad_Segments)<0);
    
    % Correct for markers at BOF & EOF:
    if length(OnsetsBadIdx)>length(OffsetsBadIdx)
        OffsetsBadIdx(end+1) = size(EEGd,2); %#ok<SAGROW>
    end
    if length(OffsetsBadIdx)>length(OnsetsBadIdx)
        OnsetsBadIdxBKP = OnsetsBadIdx;
        OnsetsBadIdx = nan(size(OnsetsBadIdx)+[1 0]);
        OnsetsBadIdx(2:end) = OnsetsBadIdxBKP;
        OnsetsBadIdx(1) = 1;
    end
    if OnsetsBadIdx(1)>OffsetsBadIdx(1)
        % In this case, both BOF and EOF are bad segments:
        OnsetsBadIdxBKP = OnsetsBadIdx;
        OnsetsBadIdx = nan(size(OnsetsBadIdx)+[1 0]);
        OnsetsBadIdx(2:end) = OnsetsBadIdxBKP;
        OnsetsBadIdx(1) = 1;
        
        OffsetsBadIdx(end+1) = size(EEGd,2); %#ok<SAGROW>
    end
    
    % Construct additional markers:
    LabBad = repmat({'"scalp_BAD"'},length(OnsetsBadIdx),1);
    % Overwrite previous .mrk file
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],...
        [MstartD-1;FileJunctionsNewFreq.start(2:end)-1;OnsetsBadIdx-1],...
        [MendD-1;FileJunctionsNewFreq.start(2:end)-1;OffsetsBadIdx-1],...
        [Mlab;repmat({'"File_junction"'},length(FileJunctionsNewFreq.start)-1,1);LabBad]);
    write_mrk_file_AnyWave([fullfile(OutFilePath,OutFilename),'.ades.mrk'],...
        [Mlab;repmat({'"File_junction"'},length(FileJunctionsNewFreq.start)-1,1);LabBad],...
        [MstartD-1;FileJunctionsNewFreq.start(2:end)-1;OnsetsBadIdx-1],...
        [MendD-1;FileJunctionsNewFreq.start(2:end)-1;OffsetsBadIdx-1]-[MstartD-1;FileJunctionsNewFreq.start(2:end)-1;OnsetsBadIdx-1]);
    
    % % for visualization in Cartool:
    % clipboard('copy',num2str(Bad_Channels'))
    % clipboard('copy',num2str(MstartDfaster'))
    
    % Save results
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_artefacts'];
    
    fprintf('================================\n================================\nArtefacts detection done, saving...\n================================\n================================\n')
    save(fullfile(OutFilePath,[OutFilename,'.mat']),...
        'Bad_Channels', 'Good_Channels', 'Bad_Segments', 'Good_Segments',...
        'EGI', 'Metrics', 'ArtDetectResults', 'OnsetsBadIdx', 'OffsetsBadIdx')
    
    % #TODO: try this pipeline on SEEG data as well...? Spikes  will likely
    % induce false positives... Anyhow, there are already some markers of nice
    % rest... Maybe check later F-tract's pipeline... (for bad channels
    % detection at least) !
    
    %% ICA
    % In contrast to FASTER pipeline, we perform ICA on whole EEG, with segments
    % of bad EEG rejected based on sliding-window approach combined with
    % original FASTER pipeline, but we do not do baseline correction, and we
    % do not interpolate electrodes before ICA, because redundant information
    % does not need to be decomposed. Instead, we ensure that ICA is performed
    % on good channels only, possibly less channels, which might be better in
    % the end, and re-include missing channels by interpolating them after.
    % Also, we won't further interpolate bad channels within epochs,
    % because either a channel is constantly bad throughout the recording,
    % or it is bad for only a specific time period (segment), but as we
    % detect these before updating bad channels, interpolating different
    % electrodes per epoch should in principle not be needed...
    fprintf('================================\n================================\nStarting ICA decomposition...\n================================\n================================\n\n')
    %=========== ICA decomposition ===========
    [W,Winv,Activations,Topos,CompWithCardiacFreq,FreqPower,FreqList] = ...
        my_ICA( EEGd(Good_Channels,Good_Segments), NewFreq );
    
    fprintf('================================\n================================\nICA decomposition done, saving...\n================================\n================================\n')
    % Save results
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_ICA'];
    
    save(fullfile(OutFilePath,[OutFilename,'.mat']),'W','Winv','Activations','Topos','CompWithCardiacFreq','FreqPower','FreqList')
    
    % Adjust markers after bad segments excision:
    [AdjMstartD,AdjMendD,AdjMlab] = adjust_mrk_4_bad_segments(MstartD,MendD,Mlab,OnsetsBadIdx,OffsetsBadIdx);
    
    % Write XYZ without bad channels for visualization in Cartool:
    write_xyz_file([OutFilePath,filesep,'Template_cap_without_bad_electrodes.xyz'],EGI.X(Good_Channels),EGI.Y(Good_Channels),EGI.Z(Good_Channels),EGI.clinicalname(Good_Channels))
    
    % Write activations time course and topography for visualization in Cartool:
    CompNames = regexprep(cellstr([repmat('C',size(Activations,1),1),num2str((1:size(Activations,1))')]),' ','');
    write_sef(fullfile(OutFilePath,[OutFilename,'_time_courses.sef']),Activations',NewFreq,CompNames);
    % Below we add zeros for the first "time frame" because Cartool starts
    % at 0. Otherwise component 1 is at time frame 0, component 2 at time
    % frame 1, ... and so on... This avoids making stupid errors when
    % looking at IC topography:
    write_sef(fullfile(OutFilePath,[OutFilename,'_topo.sef']),[zeros(1,size(Topos,2));Topos],NewFreq,ChannelNames(Good_Channels));
    write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'_time_courses.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % to directly open ICs topography with appropriate XYZ:
    write_lm_file(fullfile(OutFilePath,[OutFilename,'_topo.lm']),[{fullfile(OutFilePath,[OutFilename,'_topo.sef'])};{[OutFilePath,filesep,'Template_cap_without_bad_electrodes.xyz']}]);
    
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
    write_sef(fullfile(OutFilePath,[OutFilename,'_vs_EOG.sef']),CompWithEOG',NewFreq,CompWithEOGlab);
    write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'_vs_EOG.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    CompEMGcorr = corr(Activations',EMGlike');
    [~,CompEMGcorrMidx] = sort(abs(mean(CompEMGcorr,2)),'descend'); % could also be mean(abs()) but a given IC should in principle be correlated the same way with all bipolar channels because they share very similar info
    CompWithEMG = [EMGlike;Activations(CompEMGcorrMidx,:)]; % Good_Segments not needed here, bad segments already excised
    CompWithEMGlab = [regexprep(cellstr([char(Chan4emg(:,1)),repmat('-',size(Chan4emg,1),1),char(Chan4emg(:,2))]),' ','');CompNames(CompEMGcorrMidx)];
    write_sef(fullfile(OutFilePath,[OutFilename,'_vs_EMG.sef']),CompWithEMG',NewFreq,CompWithEMGlab);
    write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'_vs_EMG.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
    
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
    
        %========= Write ICA results that can be loaded in AnyWave =========
    
    %========= review ICA components =========
    % #TODO: pass ICAcorr as argin to review_ICA and make plot like above
    % #TODO: open not only activations and EEG, but make .ades file and open it
    % from AnyWave with these variables in ...XXX.ades.ica.mat file:
    %
    %                     hpf: [1x1     double]
    %                  labels: [1x254   cell]     => should match labels defined in .ades file for scalp !
    %                     lpf: [1x1     double]
    %                  mixing: [254x120 double]
    %                modality: [1x3     char]
    %                      sr: [1x1     double]
    %                unmixing: [120x254 double]
    
    % => AnyWave crashes, but I will use my custom function to check for the
    % effect of removing a given component! EDIT 2019-06-17 #RM@FBMlab: the
    % crash was caused by unmatched labels (see below), now it is solved !
    
        OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
        TempOutFilename = regexp(OutFilename,'_','split');
        OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_ICA_AnyWave'];
    
        hpf = HighPassBand; %%#ok<NASGU>
        lpf = LowPassBand; %%#ok<NASGU>
%         labels = ChannelNames(Good_Channels); % CRASHES ANYWAVE BECAUSE
% %         DOES NOT MATCH LABELS DEFINED IN .ADES EARLIER !
        labels = EGI.clinicalname(Good_Channels); %%#ok<NASGU>
        mixing = Winv; unmixing = W; %%#ok<NASGU>
        modality = 'EEG'; %%#ok<NASGU>
        sr = NewFreq; %%#ok<NASGU>
    
        save(fullfile(OutFilePath,[OutFilename,'.mat']),'hpf','lpf','labels','mixing','unmixing','modality','sr')
        clear hpf lpf labels mixing unmixing modality sr
    
    %======= Write .sef with only good segments and good channels for ICA review =======
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz_only_good_chan_and_TF'];
    write_sef(fullfile(OutFilePath,[OutFilename,'.sef']),EEGd(Good_Channels,Good_Segments)',NewFreq,vect(ChannelNames(Good_Channels)));
    write_mrk_file_Cartool(fullfile(OutFilePath,[OutFilename,'.sef.mrk']),AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    warning('Now waiting for user to inspect ICs !...')
    
end
diary off

%% ========== MANUAL PART ==========
%%=========== ICs review ===========

subj=1
subj=2
subj=3
subj=4
% subj=5
% subj=6

[~,~,ScalpSefFilepaths] = my_recursive_listfiles(fullfile(RootPath,SubPath{subj},DataType),'hdEEG');
ScalpSefFilepaths = cellstr(ScalpSefFilepaths);
ScalpSefFilepaths = ScalpSefFilepaths(~cellfun(@isempty,regexp(ScalpSefFilepaths,'.sef$')));
OutFilePath = regexprep(spm_file(spm_file(ScalpSefFilepaths{1},'path'),'path'),'2b_Realigned','3_Preprocessed');
OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
TempOutFilename = regexp(OutFilename,'_','split');
% Give summary of bad channels & segments detection:
OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_artefacts'];
load(fullfile(OutFilePath,[OutFilename,'.mat']))
fprintf('\n%d channels were flagged as bad...\n',length(Bad_Channels))
disp(EGI.clinicalname(Bad_Channels))
fprintf('\n%.2f %% of time frames were marked as bad...\n',sum(Bad_Segments)/length(Bad_Segments)*100)
% Load ICA results:
OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_ICA'];
load(fullfile(OutFilePath,[OutFilename,'_infos.mat']))
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~ Review ICs ~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
review_ICA(fullfile(OutFilePath,[OutFilename,'.mat']),...
    fullfile(OutFilePath,[OutFilename,'_infos.mat']),NewFreq,1:(15*60*NewFreq));


review_ICA(fullfile(OutFilePath,[OutFilename,'.mat']),...
    fullfile(OutFilePath,[OutFilename,'_infos.mat']),NewFreq);


% Not first 10 minutes but after, e.g. 2nd to 12th minute:
review_ICA(fullfile(OutFilePath,[OutFilename,'.mat']),...
    fullfile(OutFilePath,[OutFilename,'_infos.mat']),NewFreq,(2*60*NewFreq):(12*60*NewFreq));



%======= Useful: =======
clipboard('copy',num2str(Bad_Channels'))
%=======================

%======= Might be useful: =======
Good_SegmentsIdx = find(Good_Segments);
TimePointInAllSegments = Good_SegmentsIdx(TimePointInGoodSegments)
%================================

%% =========== ADDITIONAL DEFINITIONS FOR NEXT PART ============
%========== ICA components to remove for each subject ==========
IC2r = {[1,3,4,5,10,21,41,51,60,66,84,97,112] % cardiac: 10; ocular: 1,3,4; ocular/muscular: 5; non-physiological: 21,41,51,60,66,84,97,112
        [10,18,19,92] % cardiac: 18; ocular: 10,19,92
        [23,33] % cardiac: 33; ocular: 23
        [1,2,3,5,13,25,32,49,63,82,88,89,91]}; % cardiac: 13; ocular: 1,2,3,5; non-physiological: 25,32,63,88,49,89,82,91
% IC2r = {[9,13,17,19,20,24,26,32,46,48,49,50,65,70,80,81,83,93,95,100,101] % ocular: 13,9; cardiac: 19; muscular: not found; noisy electrode: 46, 49; electrode pop: 32, 70, 80, 93, 95, 100, 101; electrode drift: 17, 20, 24, 26, 48, 50, 65, 81, 83
%     [1,4,9,15,23,30,39,49,50,62,63,64] % ocular: 1, 4, 9, 15, 23, 30, 64; cardiac: 50 & 62 (artefact affected different electrodes during recording); muscular: not found; noisy electrode: not found; electrode pop: 39, 49 ; electrode drift: 63
%     [5,16,37,38,45,50,58,87,97] % ocular: 5, 16, 37, 38; cardiac: 50; muscular: 58 ; noisy electrode: not found ; electrode pop: 45, 87, 97 ; electrode drift: not found
%     [3,17,24,51] % ocular: 3 ; cardiac: 17 ; muscular: not found; noisy electrode: not found ; electrode pop: 24, 51 ; electrode drift: not found
%     [2,5,16,19,21,22,24,28,30,32,34,35,36,37,39,40,46,47,52,55,57,58,59,60,62,64,65,66,67,71,72,76,79,80,82,95,96,99,102] % ocular: 2, 5, 16, 19 ; cardiac: 32 ; muscular: ; noisy electrode: 24 ; electrode pop: 22, 35, 37, 39, 47, 55, 57, 58, 60, 62, 65, 66, 67, 71, 72, 79, 80, 82, 95, 96, 102 ; electrode drift: 28, 30, 34, 36, 40, 46, 52, 59, 64, 21, 76, 99
%     [3,12,13,18,27,35,36,40,42,47,52,53,57,66,70,81,86]}; % ocular: 3, 12 ; cardiac: 13, 27 ; muscular: ; noisy electrode: ; electrode pop: 35, 36, 40, 47, 52, 53, 57, 66, 70, 81, 86 ; electrode drift: 18, 42

XYZSubjSpacePath = {'E:\FS_subjects_DONE\sub-50\LSMAC_no_cereb\More\LSMAC_no_cereb.Original.xyz'
'E:\FS_subjects_DONE\sub-48\LSMAC_no_cereb\More\LSMAC_no_cereb.Original.xyz'
'E:\FS_subjects_DONE\sub-82\LSMAC_no_cereb\More\LSMAC_no_cereb.Original.xyz'
'E:\FS_subjects_DONE\sub-96\LSMAC_no_cereb\More\LSMAC_no_cereb.Original.xyz'};

%% ================================= PART II ===================================
%=========== Load again all files per subject and proceed to part II ===========

diary(fullfile(PreprocPath,filesep,'Log_pipeline_v1.txt'))
cowsay('Continuing EEG processing pipeline V1...')
for subj = SubjInd
    
    fprintf('%s\n',datestr(now,'yyyy.mm.dd, hh:MM:ss','local'))
    
    fprintf('Doing subject %d/%d...\n',subj,length(SubjInd))
    
    fprintf('================================\n================================\nLoading back all files from first part...\n================================\n================================\n')
    % load all files:
    [~,~,ScalpSefFilepaths] = my_recursive_listfiles(fullfile(RootPath,SubPath{subj},DataType),'hdEEG');
    ScalpSefFilepaths = cellstr(ScalpSefFilepaths);
    ScalpSefFilepaths = ScalpSefFilepaths(~cellfun(@isempty,regexp(ScalpSefFilepaths,'.sef$')));
    SEEGsefFilepaths = regexprep(ScalpSefFilepaths,'hdEEG_','icEEGi_');
    [~,~,XYZfidPath] = my_recursive_listfiles(fullfile(PreprocPath,SubPath{subj},SubPathXYZnorm),'.ToFiducial.xyz');
    
    OutFilePath = regexprep(spm_file(spm_file(ScalpSefFilepaths{1},'path'),'path'),'2b_Realigned','3_Preprocessed');
    
    % % icEEG is not needed anymore:
    %     OutFilename = spm_file(SEEGsefFilepaths{1},'basename');
    %     TempOutFilename = regexp(OutFilename,'_','split');
    %     OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    %     load([OutFilePath,filesep,OutFilename,'.mat']);
    
    OutFilename = spm_file(ScalpSefFilepaths{1},'basename');
    TempOutFilename = regexp(OutFilename,'_','split');
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    load([OutFilePath,filesep,OutFilename,'.mat']);
    
    % % ECG is not needed anymore:
    %     OutFilename = ['cardiac_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    %     load([OutFilePath,filesep,OutFilename,'.mat']);
    %     OutFilename = ['cardiac_',TempOutFilename{2},'_',TempOutFilename{3},'_bp1-15hz_down',num2str(NewFreq),'hz'];
    %     load([OutFilePath,filesep,OutFilename,'.mat']);
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_bp',num2str(HighPassBand),'-',num2str(LowPassBand),'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz'];
    load([OutFilePath,filesep,OutFilename,'.mat']);
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_artefacts'];
    load(fullfile(OutFilePath,[OutFilename,'.mat']))
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},'_ICA'];
    load(fullfile(OutFilePath,[OutFilename,'.mat']))
    load(fullfile(OutFilePath,[OutFilename,'_infos.mat']),'ICAinfos');
    
    %% ICA reconstruction (back-transform)
    
    fprintf('================================\n================================\nICA reconstruction...\n================================\n================================\n')
    ICsRemoved = IC2r{subj};
    ICAresultsFile = fullfile(OutFilePath,[OutFilename,'.mat']);
    fprintf('Components to remove for subject %d:\n', subj)
    disp(ICsRemoved')
    
    EEGc = my_ICA_backtransform(fullfile(OutFilePath,[OutFilename,'.mat']),IC2r{subj});
    
    % Re-adjust markers timings
    [AdjMstartD,AdjMendD,AdjMlab] = adjust_mrk_4_bad_segments(MstartD,MendD,Mlab,OnsetsBadIdx,OffsetsBadIdx);
    
    fprintf('================================\n================================\nSaving output...\n================================\n================================\n')
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_chan_and_TF'];
    ChannelNamesOK = ChannelNames(Good_Channels);
    NchansOK = length(ChannelNamesOK);
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGc','NewFreq','ChannelNamesOK','NchansOK',...
        'ICsRemoved','ICAresultsFile',...
        'AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGc',NewFreq,ChannelNamesOK);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(Good_Channels,Good_Segments) = EEGc;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_chan_and_TF'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'ICsRemoved','ICAresultsFile',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,ChannelNames);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    %% Bad channel interpolation
    % (after ICA, because redundant information does not need to be decomposed!)
    fprintf('================================\n================================\nBad channel interpolation using 3D spline linear interpolation...\n================================\n================================\n')
    % 3D spline (linear) interpolation
    % We will use the template average XYZ, which should be sufficient for what
    % we are aiming at here, and avoids having to center and normalize the
    % subject-specific XYZ (using e.g. normalize_center_xyz_Cartool.m)
    sXYZ = get_centered_normalized_257ch_EGI_template;
    % XYZ = sXYZ; % here we will first try with normalized centered template, but later we will extend it to subject-specific normalized & centered electrode coordinates
    
    % to get the subject-specific XYZ we anyhow need to use Cartool beforehand,
    % so we will compute it at the same time we perform ESI (instead of using
    % normalize_center_xyz_Cartool.m):
    [X,Y,Z] = read_xyz_file(XYZfidPath);
    XYZfid = [X,Y,Z];
    
    % For this function, we need the scaled / normalized / centered XYZ
    % (created by Cartool) to be passed to the 3D spline:
    %     EEGi =
    %     interp_3D_lin_bad_elec(EEGfullNaNzeroed(:,Good_Segments),XYZfid,Bad_Channels,'scattered_interpolant');
    %     % too slow...
    EEGi = interp_3D_lin_bad_elec(EEGfullNaNzeroed(:,Good_Segments),XYZfid,Bad_Channels,'polyharmonic_splines');
    
    fprintf('================================\n================================\nSaving output...\n================================\n================================\n')
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(:,Good_Segments) = EEGi;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_TF_i3DS2'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGi','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGi',NewFreq,ChannelNames);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_TF_i3DS2'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,ChannelNames);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    %% Spatial Interseptile Weighted Mean
    fprintf('================================\n================================\nApplying spatial interspeptile weighted mean filter...\n================================\n================================\n')
    
    [X,Y,Z,Names] = read_xyz_file(XYZSubjSpacePath{subj});
    XYZss = [X,Y,Z];
    % Here, on the contrary, we need the original XYZ in subject-space, because
    % the central electrodes will be given a weight of 1: if we pass the scaled
    % / normalized / centered XYZ, the central electrode will be penalized
    % because neighbours will have bigger weights! Therefore we need the XYZ
    % with the original coordinates in mm, such that inverse distance weighting
    % works as it should (NB: rotations & translations won't affect distances):
    EEGisf = spatial_interseptile_weighted_mean(EEGi,XYZss);
    EEGgsf = spatial_interseptile_weighted_mean(EEGc,XYZss(Good_Channels,:));
    
    fprintf('================================\n================================\nSaving output for spatially filtered interpolated and recontruscted EEG...\n================================\n================================\n')
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(:,Good_Segments) = EEGisf;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_TF_i3DS2_SIWMF'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGisf','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGisf',NewFreq,ChannelNames);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_TF_i3DS2_SIWMF'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,ChannelNames);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    fprintf('================================\n================================\nSaving output for spatially filtered reconstructed EEG without bad channels...\n================================\n================================\n')
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(Good_Channels,Good_Segments) = EEGgsf;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_chan_and_TF_SIWMF'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGisf','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGgsf',NewFreq,ChannelNames(Good_Channels));
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_TF_and_good_chan_SIWMF'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,ChannelNames);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    %% Re-referencing to average
    fprintf('================================\n================================\nRe-referencing to average...\n================================\n================================\n')
    % Re-reference spatially filtered interpolated (and reconstructed) EEG:
    EEGisfr = bsxfun(@minus,EEGisf,mean(EEGisf,1));
    % ... and re-reference also EEG with spatial filter but no interpolation:
    EEGgsfr = bsxfun(@minus,EEGgsf,mean(EEGgsf,1));
    % ... and re-reference also EEG without spatial filter but interpolated:
    EEGir = bsxfun(@minus,EEGi,mean(EEGi,1));
    % ... and re-reference also EEG without spatial filter and interpolation:
    EEGgr = bsxfun(@minus,EEGc,mean(EEGc,1));
    
    fprintf('================================\n================================\nSaving output (re-referenced spatially filtered interpolated EEG)...\n================================\n================================\n')
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(:,Good_Segments) = EEGisfr;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_TF_i3DS2_SIWMF_avgref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGisfr','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGisfr',NewFreq,EGI.clinicalname);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_TF_i3DS2_SIWMF_avg_ref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,EGI.clinicalname);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    fprintf('================================\n================================\nSaving output (re-referenced EEG spatially filtered without bad channels)...\n================================\n================================\n')
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(Good_Channels,Good_Segments) = EEGgsfr;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_chan_and_TF_SIWMF_avgref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGgsfr','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGgsfr',NewFreq,EGI.clinicalname(Good_Channels));
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_TF_but_good_chan_SIWMF_avg_ref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,EGI.clinicalname);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    fprintf('================================\n================================\nSaving output (re-referenced interpolated EEG)...\n================================\n================================\n')
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(:,Good_Segments) = EEGir;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_TF_i3DS2_avgref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGir','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGir',NewFreq,EGI.clinicalname);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_TF_i3DS2_avg_ref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,EGI.clinicalname);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    fprintf('================================\n================================\nSaving output (re-referenced EEG without bad channels)...\n================================\n================================\n')
    % re-allocate inside whole:
    EEGfull = nan(size(EEGd)); EEGfull(Good_Channels,Good_Segments) = EEGgr;
    EEGfullNaNzeroed = EEGfull; EEGfullNaNzeroed(isnan(EEGfullNaNzeroed))=0;
    
    % Write output without bad segments:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_only_good_chan_and_TF_avgref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGgr','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','AdjMstartD','AdjMendD','AdjMlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGgr',NewFreq,EGI.clinicalname(Good_Channels));
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],AdjMstartD-1,AdjMendD-1,AdjMlab);
    
    % Write output with bad segments NaNed / zeroed:
    OutFilename = [TempOutFilename{1},'_',TempOutFilename{2},'_',TempOutFilename{3},...
        '_bp',num2str(HighPassBand),'-',num2str(LowPassBand),...
        'hz_notch',num2str(LineNoise),'_down',num2str(NewFreq),'hz','_ICA_recon_with_bad_TF_but_good_chan_avg_ref'];
    save([OutFilePath,filesep,OutFilename,'.mat'],...
        'EEGfull','NewFreq','ChannelNames','Nchans',...
        'Bad_Channels','Bad_Segments','OnsetsBadIdx','OffsetsBadIdx',...
        'MstartD','MendD','Mlab','-v7.3');
    write_sef([OutFilePath,filesep,OutFilename,'.sef'],...
        EEGfullNaNzeroed',NewFreq,EGI.clinicalname);
    write_mrk_file_Cartool([OutFilePath,filesep,OutFilename,'.sef.mrk'],MstartD-1,MendD-1,Mlab);
    
    %% Make epochs, concatenate, run analyses... #TODO!
    
    % => outside: need .mat files!
    
    
    % ...remember to adjust marker timings based on Good_Segments! => not
    % needed anymore, we adjusted the timings before downsampling!
    
    %% #TODO: perform microstates analysis in Matlab (EEGLAB? Brainstorm?) as well?...
    
end
diary off

%% =========== CHECK MARKERS FOR DICTIONARIES ============

%========== Check that all markers of a given type have duration = 0 or > 0 ==========
% Construct dictionaries:
subj=1 % #OK
subj=2 % #OK
subj=3 % #TOFIX
subj=4 % #OK
subj=5 % #TOFIX
subj=6 % #TOFIX
% PreprocFile = my_listfiles(fullfile(PreprocPath,SubPath{subj},DataType),'with_bad_TF_i3DS2_SIWMF_avg_ref.mat');
% PreprocFile = char(fullfile(PreprocPath,SubPath{subj},DataType,PreprocFile));
PreprocFile = spm_select('FPList',fullfile(PreprocPath,SubPath{subj},DataType),'^hdEEG_.*down200hz\.mat$');
Table = dual_check_mrk(cellstr(spm_file(PreprocFile,'ext','sef.mrk')),1);
for mrktype = 1:size(Table,1)
    % There is only 1 file each time, so we can hardcode cell in cell to 1:
    DurTest = Table{mrktype,4}{1}-Table{mrktype,3}{1};
    if min(DurTest)==0 && max(DurTest)~=0
        warning('Possible issue with marker %s!',Table{mrktype,1})
        fprintf('\n%s:\nn=%d, min=%d, avg=%d, max=%d\n',Table{mrktype,1},str2double(Table{mrktype,2}),min(DurTest),mean(DurTest),max(DurTest))
        fprintf('Number of durations > 0: %d\n',sum(DurTest>0))
        disp(Table{mrktype,3}{1}(DurTest>0))
    end
end
Marker_Times = [173472
    341484];
ChannelOfInterest = {'HPD1-HPD2'};
SEEGfile = spm_select('FPList',fullfile(PreprocPath,SubPath{subj},DataType),'^icEEGi_.*down200hz\.mat$');
load(SEEGfile,'icEEGd','icChannelNames')
[~,~,labelsB,~,EEGtemp] = bipolar_montage(icChannelNames,2,icEEGd);
ChannelToShow = match_vectors(ChannelOfInterest,labelsB,1);
[Shifts,ToCheck] = spike_aligner_v3( EEGtemp, Marker_Times, ChannelToShow, ChannelOfInterest, 200 );
% Mrk = Mrk-Shifts % that works, has been tested.

% When I previously used spike_aligner_wrapper_v2.m to align spikes,
% in spike_aligner_icEEG.m, I did this stupid thing of using the average
% between marker start and marker end, like this:
% Marker_Time = round((Marker_T1(MarkersToDo)+Marker_T2(MarkersToDo))/2);
% That means that the average between the two time frames was manually
% aligned. Thus, if there are still markers with durations > 0 and that
% should have duration = 0, I can take again the rounded average and
% consider it as being aligned. Therefore, there is no need (at least in
% the present case (i.e. for these 6 patients), to check again alignment to
% determine whether what is the time frame that is aligned. I just take the
% rounded average again. And because I did this for all important marker
% types (i.e. all spikes markers), even e.g. the "mini" are aligned. I just
% need to figure out which are spike markers and which are not, but for
% this, I can just look again at spike_aligner_wrapper_v2.m to build the
% dictionary and determine which are the "less important" spikes which I
% should consider when purifying epochs.

%% ========= To check dictionary against markers =========
subj=1
subj=2
subj=3
subj=4
subj=5
subj=6
PreprocFile = my_listfiles(fullfile(PreprocPath,SubPath{subj},DataType),'with_bad_TF_i3DS2_SIWMF_avg_ref.mat');
PreprocFile = char(fullfile(PreprocPath,SubPath{subj},DataType,PreprocFile));
Table = dual_check_mrk(cellstr(spm_file(PreprocFile,'ext','sef.mrk')),1) % visualize counts & get more infos
dual_check_mrk(cellstr(spm_file(PreprocFile,'ext','sef.mrk')),1) % visualize counts only
Table = dual_check_mrk(cellstr(spm_file(PreprocFile,'ext','sef.mrk')),0) % get more infos only

% when checking multiple subjects, to know which is which...
title(SubPath{subj})

% To check whether some markers co-occur:
Closest = match2bins(Table{5,3}{1},Table{8,3}{1});
% figure; scatter(Table{5,3}{1},Table{8,3}{1}(Closest))
Table{5,3}{1}-Table{8,3}{1}(Closest)

%% ==================== To build exceptions ====================
% To check for spike markers of the same family:
spk=1
spk=2
spk=3
spk=4
spk=5
spk=6
spk=7
ThisSpike = regexprep(SpikeLabel{subj}{spk},'\^|\$|"','');
ThisSpike = regexp(ThisSpike,'_','split');

ThisSpike = ThisSpike{1};
ThisSpike = ThisSpike{2}; % for RD...

FamilyTemp = regexprep(Table(~cellfun(@isempty,regexp(Table(:,1),ThisSpike)),1),'^("SEEG_|"scalp_)','"');
for spkfam = 1:length(FamilyTemp)
    if strcmp(FamilyTemp{spkfam}(1),'"') && strcmp(FamilyTemp{spkfam}(end),'"')
        FamilyTemp{spkfam} = ['^"',FamilyTemp{spkfam}(2:end-1),'"$'];
    else
        error('Some spike markers do not start / end with ''"'' !')
    end
end

%% =========== ADDITIONAL DEFINITIONS FOR NEXT PART ============
% Adapted from spike_aligner_wrapper_v2.m:
% In what follows, we ignored "Epitome" markers, "?" and other non-spike
% markers. We ignore mini markers and non-spikes events such
% as rhythmic activity, but we avoid other spike types
% (those not in SEEGmarkersList won't be considered for exclusion at all
% when looking at a specific case in SpikeLabel):

% NB: the second column, which contains electrode labels involved for each
% spike, is actually omre useful for spike alignment (refinement) and for
% future analyses but not for epoching per se).

% "ARTEFACT" markers should not be retained here, because we only care
% about scalp EEG quality, which has been already assessed earlier (see
% part I of pipeline)

%======== icEEG ========
% AY:
SEEGmarkersList{1,1} = {'^"CP1"$|^"CP2"$' % occurs systematically 200-300 ms after HAD/HPD
                        '^"HAD"$|^"HAD2"$' % only HAD
                        '^"PT3"$'
                        '^"PT6"$'
                        '^"bHAD_HPD"$|^"bHPD_HAD"$|^"HAD_HPD"$'
                        '^"mHAD"$'
                        '^"mHPD"$'};
% '^"bHAD_HPD"$|^"bHPD_HAD"$': merging all conditions, consider for averaging
% '^"HAD"$|^"HAD2"$': consider for averaging
% Other spikes are too small or rare to consider them...

% BY
SEEGmarkersList{2,1} = {'^"AG11"$|^"ag11"$' % lowercase vs uppercase...
                        '^"ETD"$|^"ETD9"$' % only the one with number mentioned in notes, guess the other relates to the same event
                        '^"ETG"$|^"ETG5"$' % only the one with number mentioned in notes, guess the other relates to the same event
                        '^"HAD9"$'
                        '^"HAD_HPD"$'
                        '^"bHAG_HPG"$|^"HAG_HPG"$|^"HAG"$'
                        '^"mETG_HPG_HAG"$'
                        '^"mHAG_HPG"$'
                        '^"nbHAG"$'};
% '^"bHAG_HPG"$|^"HAG_HPG"$|^"HAG"$': merging all mentioned categories,
% consider for averaging
% '^"ETD"$|^"ETD9"$': consider for averaging
% '^"ETG"$|^"ETG5"$': consider for averaging
% Other spikes are too small or rare to consider them...

% CC
SEEGmarkersList{3,1} = {'^"AD"$' % not mentioned by SV in notes, rare
                        '^"AG"$' % from SV's notes, not clear whether it can be used yet or not...
                        '^"AG_HAG"$|^"AG_HAG_HPG"$|^"AG_HAG_HPGlat"$|^"AG_HAG_PTG"$' % second not mentioned in the notes, although first mention lateral contacts, so likely the same, third merged as suggested by SV in notes
                        '^"FAD"$' % infrequent, small, medial
                        '^"FAD9"$' % infrequent, small, lateral
                        '^"HAG"$' % rare
                        '^"HAGlat_HPGlat"$|^"HPGlat"$|^"HPG"$' % the two last are not mentioned in the notes by SV, but a quick look at the averages suggests I could consider them as being the same event...
                        % '^"HPD_HFO"$'
                        '^"IDG"$|^"IDG5"$' % infrequent, small
                        '^"IVG"$|^"IVGlat"$' % infrequent, small
                        '^"TGlat"$' % rare
                        '^"bHAD_HPD"$|^"bbHAD_HPD_AD"$|^"HAD_HPD_AD"$|^"HAD_HPG_AD"$|^"HAD_HPD"$' % "bbHAD_HPD" & "bHAD_HPD" mentioned in the notes by SV, but "AD" part in "bbHAD_HPD" not explicit, I guess I can merge them because there are few of each category... "HAD_HPG_AD" not mentioned in notes, likely mistyped by SV...
                        '^"burstFAD"$' % rare
                        '^"burstIVGlat"$' % rare
                        '^"mAG"$' % small
                        '^"mAG_HAG"$' % rare and small
                        '^"mAG_HAG_HPG"$' % small
                        '^"mAGneg"$' % rare and small
                        '^"mAGpos"$' % rare and small
                        '^"mHAD_HPD"$' % frequent but too small and too blunt
                        '^"mIDG"$' % small
                        '^"pHAD_HPD"$'}; % rare and small
% '^"AG_HAG"$|^"AG_HAG_HPG"$|^"AG_HAG_HPGlat"$|^"AG_HAG_PTG"$': consider for averaging
% '^"HAGlat_HPGlat"$|^"HPGlat"$|^"HPG"$': consider for averaging...?
% '^"bHAD_HPD"$|^"bbHAD_HPD_AD"$|^"HAD_HPD_AD"$|^"HAD_HPG_AD"$|^"HAD_HPD"$': consider for averaging! 
% Other spikes are too small or rare to consider them...

% MP
SEEGmarkersList{4,1} = {'^"CA23rhythmictheta"$|^"CAD23"$|^"CADrhzthmictheta"$|^"CAD_IAD"$' % 4th not mentioned by SV in notes, I guess it relates to the same kind of event with more IAD component... Looking at averages, that's what it suggests...
                        '^"HAD_HPD"$'
                        '^"HAD_HPD_AD"$'
                        '^"HAD_HPD_HFO"$'
                        '^"bHAD_HPD"$|^"bbHAD_HPD_AD"$|^"bHAD_HPD_AD"$' % merged based on similarity of averages
                        '^"mHAD_HPD"$'};
% '^"CA23rhythmictheta"$|^"CAD23"$|^"CADrhzthmictheta"$|^"CAD_IAD"$': consider for averaging
% '^"bHAD_HPD"$|^"bbHAD_HPD_AD"$|^"bHAD_HPD_AD"$': consider for averaging

%======== hdEEG ========
ScalpMarkersList{1,1} = {%'"scalp_qualité diminue"'
                         %'"scalp_Marker"'
                         ''};

ScalpMarkersList{2,1} = {'^"LRDA_RT"$'
                         '^"RF"$'
                         '^"TD?"$'
                         %'^"no IED, a few doubtful minithings to check on icEEG"$'
                         '^"x"$'};

ScalpMarkersList{3,1} = {'^"T7"$'
                         '^"T7map"$'
                         '^"T8"$'
                         '^"T8map"$'
                         '^"mT7"$'
                         '^"sT7"$'
                         '^"sT7map"$'};

ScalpMarkersList{4,1} = {''};

%% Definitions for spikes of interest:

SpikeLabel = {{'^"bHAD_HPD"$|^"bHPD_HAD"$|^"HAD_HPD"$' % HAD/HPD
               '^"HAD"$|^"HAD2"$'} % HAD only
    
              {'^"bHAG_HPG"$|^"HAG_HPG"$|^"HAG"$' % HAG/HPG
               '^"ETD"$|^"ETD9"$' % ETD
               '^"ETG"$|^"ETG5"$'} % ETG
    
              {'^"AG_HAG"$|^"AG_HAG_HPG"$|^"AG_HAG_HPGlat"$|^"AG_HAG_PTG"$' % AG(HAG/HPG)
               '^"HAGlat_HPGlat"$|^"HPGlat"$|^"HPG"$' % HAG/HPG
               '^"bHAD_HPD"$|^"bbHAD_HPD_AD"$|^"HAD_HPD_AD"$|^"HAD_HPG_AD"$|^"HAD_HPD"$'} % HAD/HPD
    
              {'^"CA23rhythmictheta"$|^"CAD23"$|^"CADrhzthmictheta"$|^"CAD_IAD"$' % CAD
               '^"bHAD_HPD"$|^"bbHAD_HPD_AD"$|^"bHAD_HPD_AD"$'}}; % HAD/HPD

% Patient-specific exceptions based on how the marking was done (there were
% redundancies in some cases):
Exceptions = {{{'^"CP1"$|^"CP2"$'
                '^"mHAD"$'
                '^"mHPD"$'}
                {'^"mHAD"$'} % ok for small HAD but not for posterior component, and CAD occurs preferentially with HAD/HPD so it should be avoided
                } % systematically occurs 200-300 ms after HAD/HPD
              
              {{'^"mHAG_HPG"$'
                '^"nbHAG"$'} % but not '^"mETG_HPG_HAG"$' because it might interfere with ETG/ETG5
              {''}
              {''}}
              
              {{'^"AG"$'
                '^"HAG"$'
                '^"mAG"$'
                 '^"mAG_HAG"$'
                 '^"mAG_HAG_HPG"$'
                 '^"mAGneg"$'
                 '^"mAGpos"$'}
              {'^"HAG"$'}
              {'^"mHAD_HPD"$'
              '^"pHAD_HPD"$'}}
              
              {{''}
              {'^"HAD_HPD"$'
               '^"HAD_HPD_AD"$'
               '^"HAD_HPD_HFO"$'
               '^"mHAD_HPD"$'}}
              };
% => "mini" version of "nice" spikes (i.e. chosen for analysis) should be kept
% if they appear within the epoch of a given spike, but if other spikes (and
% not only the spikes also considered for analysis) appear in the same
% window, we should exclude the epoch because it will be contaminated!
% => likewise, the polarity of the spike will matter, so if for a given
% spike, the polarity is positive, if a spike with negative polarity
% appears we should reject the epoch, because it is likely another spike
% type / different generator / source...

RandomLabel = {'^"random"$'};

%======= Parameters: epoch duration & maximal delay allowed between the two modalities =======
% Load info for initial sampling frequency (such that parameters
% defined above make sense given the downsampled signals):
OrigSampFreqFile = spm_select('FPList',fullfile(PreprocPath,SubPath{1},DataType),'^hdEEG_.*_bp1-45hz_notch50\.mat$');
if isempty(OrigSampFreqFile)
    error('Could not find file with original sampling frequency !')
else
    load(OrigSampFreqFile,'SamplingFreq')
end
% Here we assume that downsampling factor was the same for all patients
% (e.g. all were at 1000 Hz and were downsampled at 200 Hz), so only 1 file
% is necessary...
DownSamplingFactor = round(SamplingFreq/NewFreq);
ScalpSEEGmaxtimediff = fix(200/DownSamplingFactor); % 200 ms duration for slow waves, thus should be longer than that
PreStim = 999;
PostStim = 1000;
EpochIdx = unique(fix((-PreStim:PostStim)/DownSamplingFactor)); %-249:250;
% (1 second before and 1 second after, because this time we have downsampled
% at 250 Hz !)
EpochDuration = length(EpochIdx);

%% ================================= PART III ===================================

%% Epoching and concatenation (NEW)

cowsay('Epoching part of pipeline V1...')
for subj = SubjInd
    fprintf('Doing subject %d/%d...\n',subj,length(SubjInd))
    
    % Get .mat file with bad segments (better than .sef file because in
    % .mat bad segments are NaN, whereas they are zeroed in .sef):
    PreprocFile = my_listfiles(fullfile(PreprocPath,SubPath{subj},DataType),'with_bad_TF_i3DS2_SIWMF_avg_ref.mat');
    PreprocFile = char(fullfile(PreprocPath,SubPath{subj},DataType,PreprocFile));
    load(PreprocFile)
    mkdir([spm_file(PreprocFile,'path'),filesep,'Epochs'])
    
    for spk = 1:length(SpikeLabel{subj})
        
        % Get other spike markers and exclude epochs containing other events:
        % We use "SEEGmarkersList" and "ScalpMarkersList", but we also
        % keep manually marked ARTEFACTS (most scalp ARTEFACTS were
        % automatically rejected but there were some SEEG ARTEFACTS as well,
        % even if those were generally not close to other events (neither
        % spikes, nor randoms)...
        SEEG2exclude4thisSpk = match_vectors(SpikeLabel{subj}(spk),...
            SEEGmarkersList{subj,1},1); % matching current spike
        if iscell(SEEG2exclude4thisSpk)
            error('Could not match current spike to spikes list !')
        end
        SEEG2exclude4thisSpk = setdiff(SEEGmarkersList{subj,1},...
            SEEGmarkersList{subj,1}(SEEG2exclude4thisSpk)); % everything except current spike
        
        % Patient-specific exceptions:
        SEEG2exclude4thisSpk = setdiff(SEEG2exclude4thisSpk,Exceptions{subj}{spk});
        
        % Filter out "scalp" / "SEEG" at beginning of marker label
        mL = regexprep(Mlab,'^("SEEG_|"scalp_)','"');
        mT1 = MstartD;
        mT2 = MendD;
        
        try
            % SHOULD NOT CONTAIN SCALP SPIKES AT ALL, WITHIN THE ENTIRE EPOCH
            % WINDOW (otherwise it will induce a bias in the next analyses)
            [ SpkStart, SpkEnd, SpkLabel ] = purify_epochs( mT1, mT2, mL, SpikeLabel{subj}(spk), ScalpMarkersList{subj,1}, round(EpochDuration/2), 'both' );
            
            % Should not contain another SEEG spike of the other condition before
            % or after within the epoch (or spike that were not considered
            % but were marked!):
            [ SpkStart, SpkEnd ] = purify_epochs( SpkStart, SpkEnd, SpkLabel, SpikeLabel{subj}(spk), SEEG2exclude4thisSpk, round(EpochDuration/2), 'both' );
        catch
            SpkStart = [];
            SpkEnd = [];
        end
        
        %% get random timings & purify epochs:
        % In principle no random should be discarded here, because they were
        % constructed such that they do not overlap with other markers (in
        % practice there is slight overlap in very few cases, though)...
        % We won't use dual_mrk_rand_xor.m because it results in
        % variable duration epochs, so we will simply randomly select
        % epochs later, and for now we will just get all of them but this
        % time we need to filter out SpikeLabel !
        
        % SHOULD NOT CONTAIN SCALP SPIKES NEITHER, WITHIN THE ENTIRE EPOCH
        % WINDOW (otherwise it will induce a bias in the next analyses)
        try
            [ RandomStart, RandomEnd, RandomL ] = purify_epochs( mT1, mT2, mL, RandomLabel, ScalpMarkersList{subj,1}, round(EpochDuration/2), 'both' );
            
            % Below we want to exclude really all spike markers found:
            [ RandomStart, RandomEnd ] = purify_epochs( RandomStart, RandomEnd, RandomL, RandomLabel, SEEGmarkersList{subj,1}, round(EpochDuration/2), 'both' );
        catch
            RandomStart = [];
            RandomEnd = [];
        end
        
        %% Check duration and take rounded average between start and end if necessary:
        % At this stage, markers should be specific enough to check for duration:
        % In case there are markers with duration > 0, take the middle
        % timing (rounded average). This is ok based on the fact that
        % spike_aligner_wrapper_v2.m called spike_aligner_v2.m, with which
        % all spikes markers were aligend based on the average between
        % start and end:
        if any((SpkEnd-SpkStart)>0)
            warning('Some spike markers with duration > 0 have been found, this should not be the case but here it is fine to take the middle time point (rounded average) between start and end (see notes in this script about spike_aligner_wrapper_v2.m and spike_aligner_v2.m) !')
            SpkMiddle = round((SpkStart+SpkEnd)/2);
            SpkStart = SpkMiddle;
            SpkEnd = SpkMiddle;
        end
        if any((RandomEnd-RandomStart)>0)
            warning('Some random markers with duration > 0 have been found, this should not be the case but here it is fine to take the middle time point (rounded average) between start and end (see notes in this script about spike_aligner_wrapper_v2.m and spike_aligner_v2.m) !')
            RandomMiddle = round((RandomStart+RandomEnd)/2);
            RandomStart = RandomMiddle;
            RandomEnd = RandomMiddle;
        end
        
        %% Filter out markers that are too close to BOF or EOF:
        Nsamples = size(EEGfull,2);
        
        if ~isempty(RandomStart)
            % we look at start and end, just to be 100% sure
            ToDiscard = (((RandomEnd+EpochIdx(end))>Nsamples) + ...
                ((RandomEnd+EpochIdx(1))<0) + ...
                ((RandomStart+EpochIdx(end))>Nsamples) + ...
                ((RandomStart+EpochIdx(1))<0))>0;
            
            RandomStart(ToDiscard)=[];
            RandomEnd(ToDiscard)=[];
        end
        
        % we do it separately for randoms and spikes
        if ~isempty(SpkStart)
            ToDiscard = (((SpkEnd+EpochIdx(end))>Nsamples) + ...
                ((SpkEnd+EpochIdx(1))<0) + ...
                ((SpkStart+EpochIdx(end))>Nsamples) + ...
                ((SpkStart+EpochIdx(1))<0))>0;
            
            SpkStart(ToDiscard) = [];
            SpkEnd(ToDiscard) = [];
        end
        
        %% reject epochs that are too close in time (otherwise it duplicates data!):
        if ~isempty(RandomStart)
            RandomStart = prune_timings( RandomStart, EpochDuration );
            RandomEnd = prune_timings( RandomEnd, EpochDuration );
        end
        if ~isempty(SpkStart)
            SpkStart = prune_timings( SpkStart, EpochDuration );
            SpkEnd = prune_timings( SpkEnd, EpochDuration );
        end
        
        %% reject also epochs of different types that are too close in time...
        if ~isempty(SpkStart) && ~isempty(RandomStart)
            AllStart = sort([RandomStart;SpkStart]);
            AllEnd = sort([RandomEnd;SpkEnd]);
            
            AllStart = prune_timings( AllStart, EpochDuration );
            AllEnd = prune_timings( AllEnd, EpochDuration );
            
            Temp1 = match_vectors(RandomStart,AllStart,1);
            if isa(Temp1,'cell')
                RandomStart = RandomStart(~cellfun(@isempty,Temp1));
            end
            Temp2 = match_vectors(SpkStart,AllStart,1);
            if isa(Temp2,'cell')
                SpkStart = SpkStart(~cellfun(@isempty,Temp2));
            end
            Temp3 = match_vectors(RandomEnd,AllEnd,1);
            if isa(Temp3,'cell')
                RandomEnd = RandomEnd(~cellfun(@isempty,Temp3));
            end
            Temp4 = match_vectors(SpkEnd,AllEnd,1);
            if isa(Temp4,'cell')
                SpkEnd = SpkEnd(~cellfun(@isempty,Temp4));
            end
            
        end
        
        % Check if epoch contains NaN (in EEGfull, bad segments are NaNed, this
        % ensure continuous epochs !) (DO THIS BEFORE RANDOMLY SELECTING EPOCHS!):
        if ~isempty(SpkStart) && ~isempty(RandomStart)
            DiscardSpkEp = false(length(SpkStart),1);
            for ep = 1:length(SpkStart)
                EEGtemp = EEGfull(:,SpkStart(ep)+EpochIdx);
                DiscardSpkEp(ep) = any(isnan(EEGtemp(:)));
            end
            DiscardRandEp = false(length(RandomStart),1);
            for ep = 1:length(RandomStart)
                EEGtemp = EEGfull(:,RandomStart(ep)+EpochIdx);
                DiscardRandEp(ep) = any(isnan(EEGtemp(:)));
            end
            SpkStart = SpkStart(~DiscardSpkEp);
            SpkEnd = SpkEnd(~DiscardSpkEp);
            RandomStart = RandomStart(~DiscardRandEp);
            RandomEnd = RandomEnd(~DiscardRandEp);
        end
        
        %% Randomly select epochs to get same number of epochs in both conditions:
        if ~isempty(SpkStart) && ~isempty(RandomStart)
            if length(RandomStart)>length(SpkStart)
                warning('Randomly selecting epochs! All output might change!')
                RandomlySelectedRandoms = sort(randsample(1:length(RandomStart),length(SpkStart)));
                RandomStart = RandomStart(RandomlySelectedRandoms);
                RandomEnd = RandomEnd(RandomlySelectedRandoms);
            elseif length(SpkStart)>length(RandomStart)
                warning('Randomly selecting epochs! All output might change!')
                RandomlySelectedSpikes = sort(randsample(1:length(SpkStart),length(RandomStart)));
                SpkStart = SpkStart(RandomlySelectedSpikes);
                SpkEnd = SpkEnd(RandomlySelectedSpikes);
            end
            
            % make concatenated epochs:
            SpikesT = markers2epochs(SpkStart,SpkEnd,EpochIdx,'center');
            RandomsT = markers2epochs(RandomStart,RandomEnd,EpochIdx,'center');
            
            % Save epoch informations /!\ !
            save(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'suffix',['_epochs_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''),'_spike_epochs']),'SpkStart','SpikesT','SpkEnd')
            save(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'suffix',['_random_epochs_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''),'_random_epochs']),'RandomStart','RandomsT','RandomEnd')
            
            % for unconcatenated epochs:
            SpikesEpochsConcat = EEGfull(:,SpikesT);
            RandomsEpochsConcat = EEGfull(:,RandomsT);
            SpikesEpochs = reshape(SpikesEpochsConcat,size(SpikesEpochsConcat,1),EpochDuration,length(SpkStart));
            RandomsEpochs = reshape(RandomsEpochsConcat,size(RandomsEpochsConcat,1),EpochDuration,length(RandomStart));
            
            % Write .sef for segmentation in Cartool:
            write_sef(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'ext','sef','suffix',['_epochs_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                SpikesEpochsConcat',NewFreq,ChannelNames);
            
            write_sef(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'ext','sef','suffix',['_random_epochs_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                RandomsEpochsConcat',NewFreq,ChannelNames);
            
            for ep = 1:size(SpikesEpochs,3)
                write_sef(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'ext','sef','suffix',['_epoch_',num2str(ep,'%03.f'),'_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                squeeze(SpikesEpochs(:,:,ep))',NewFreq,ChannelNames);
            end
            
            for ep = 1:size(RandomsEpochs,3)
                write_sef(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'ext','sef','suffix',['_random_epoch_',num2str(ep,'%03.f'),'_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                squeeze(RandomsEpochs(:,:,ep))',NewFreq,ChannelNames);
            end
            
            % Write also the average(s):
            write_sef(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'ext','sef','suffix',['_avg_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                mean(SpikesEpochs,3)',NewFreq,ChannelNames);
            write_sef(spm_file([spm_file(PreprocFile,'path'),filesep,'Epochs',filesep,spm_file(PreprocFile,'filename')],'ext','sef','suffix',['_avg_random_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                mean(RandomsEpochs,3)',NewFreq,ChannelNames);
            
            % ... and write also average SEEG epoch...
            SEEG = spm_select('FPList',fullfile(PreprocPath,SubPath{subj},DataType),'^icEEGi_.*down.*\.mat$');
            load(SEEG,'icEEGd','icChannelNames')
            % Make bipolar montage:
            [~,~,icChannelsBip,~,SEEGbip] = bipolar_montage(icChannelNames,2,icEEGd);
            % Get SEEG TFs:
            SpikesEpochsConcatSEEG = SEEGbip(:,SpikesT);
            RandomsEpochsConcatSEEG = SEEGbip(:,RandomsT);
            SpikesEpochsSEEG = reshape(SpikesEpochsConcatSEEG,size(SpikesEpochsConcatSEEG,1),EpochDuration,length(SpkStart));
            RandomsEpochsSEEG = reshape(RandomsEpochsConcatSEEG,size(RandomsEpochsConcatSEEG,1),EpochDuration,length(RandomStart));
            
            % Write .sef for segmentation in Cartool:
            write_sef(spm_file([spm_file(SEEG,'path'),filesep,'Epochs',filesep,spm_file(SEEG,'filename')],'ext','sef','suffix',['_epochs_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                SpikesEpochsConcatSEEG',NewFreq,icChannelsBip);
            
            write_sef(spm_file([spm_file(SEEG,'path'),filesep,'Epochs',filesep,spm_file(SEEG,'filename')],'ext','sef','suffix',['_random_epochs_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                RandomsEpochsConcatSEEG',NewFreq,icChannelsBip);
            
            for ep = 1:size(SpikesEpochs,3)
                write_sef(spm_file([spm_file(SEEG,'path'),filesep,'Epochs',filesep,spm_file(SEEG,'filename')],'ext','sef','suffix',['_epoch_',num2str(ep,'%03.f'),'_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                squeeze(SpikesEpochsSEEG(:,:,ep))',NewFreq,icChannelsBip);
            end
            
            for ep = 1:size(RandomsEpochs,3)
                write_sef(spm_file([spm_file(SEEG,'path'),filesep,'Epochs',filesep,spm_file(SEEG,'filename')],'ext','sef','suffix',['_random_epoch_',num2str(ep,'%03.f'),'_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                squeeze(RandomsEpochsSEEG(:,:,ep))',NewFreq,icChannelsBip);
            end
            
            % Write average:
            write_sef(spm_file([spm_file(SEEG,'path'),filesep,'Epochs',filesep,spm_file(SEEG,'filename')],'ext','sef','suffix',['_avg_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                mean(SpikesEpochsSEEG,3)',NewFreq,icChannelsBip);
            write_sef(spm_file([spm_file(SEEG,'path'),filesep,'Epochs',filesep,spm_file(SEEG,'filename')],'ext','sef','suffix',['_avg_random_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')]),...
                mean(RandomsEpochsSEEG,3)',NewFreq,icChannelsBip);
            
            % #TODO: think of applying rm_artefacts (and ICA?) on SEEG data
            % (but anyhow neither interpolation, nor spatial interseptile
            % weighted mean filter!)
            
        else
            warning('No epochs remaining for\npatient %s, spike type %s !\n',SubPath{subj},regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''))
        end
    end
end

%% ========== Source estimation of epochs ==========
% to reconstruct scalp sources (SPs based on LSMAC):
SPIfiles = {'E:\FS_subjects_DONE\sub-50\LSMAC_no_cereb\More\T1_fixed.spi'
            'E:\FS_subjects_DONE\sub-48\LSMAC_no_cereb\More\T1_fixed.spi'
            'E:\FS_subjects_DONE\sub-82\LSMAC_no_cereb\More\T1_fixed.spi'
            'E:\FS_subjects_DONE\sub-96\LSMAC_no_cereb\More\T1_fixed.spi'};
% to reconstruct scalp sources (LAURA inverse solutions):
% NB: change "Laura" to "Loreta" to switch inverse...
ISfiles = {'E:\FS_subjects_DONE\sub-50\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
           'E:\FS_subjects_DONE\sub-48\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
           'E:\FS_subjects_DONE\sub-96\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
           'E:\FS_subjects_DONE\sub-82\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'};

% alternative including cerebellum:
SPIfilesAlt = {'E:\FS_subjects_DONE\sub-50\LSMAC_cereb\More\T1_fixed.spi'
            'E:\FS_subjects_DONE\sub-48\LSMAC_cereb\More\T1_fixed.spi'
            'E:\FS_subjects_DONE\sub-82\LSMAC_cereb\More\T1_fixed.spi'
            'E:\FS_subjects_DONE\sub-96\LSMAC_cereb\More\T1_fixed.spi'};
ISfilesAlt = {'E:\FS_subjects_DONE\sub-48\LSMAC_cereb\LSMAC_cereb.Laura.is'
           'E:\FS_subjects_DONE\sub-96\LSMAC_cereb\LSMAC_cereb.Laura.is'
           'E:\FS_subjects_DONE\sub-82\LSMAC_cereb\LSMAC_cereb.Laura.is'
           'E:\FS_subjects_DONE\sub-50\LSMAC_cereb\LSMAC_cereb.Laura.is'};

% for SEEG:
ELSfiles = {'E:\FS_subjects_DONE\sub-50\elec_recon\patient50_align_equ.els'
            'E:\FS_subjects_DONE\sub-48\elec_recon\patient48_align_equ.els'
            'E:\FS_subjects_DONE\sub-82\elec_recon\patient82_align_equ.els'
            'E:\FS_subjects_DONE\sub-96\elec_recon\patient96_align_equ.els'};


%% ========== Count number of epochs per spike and patient ==========
clear Nepochs
Count = 0;
for subj = SubjInd
    PreprocFile = my_listfiles(fullfile(PreprocPath,SubPath{subj},DataType),'with_bad_TF_i3DS2_SIWMF_avg_ref.mat');
    PreprocFile = char(fullfile(PreprocPath,SubPath{subj},DataType,PreprocFile));
    for spk = 1:length(SpikeLabel{subj})
        Count = Count+1;
        SpkFiles = cellstr(spm_select('FPList',[fileparts(PreprocFile),filesep,'Epochs'],spm_file(spm_file(PreprocFile,'basename'),'ext','sef','suffix',['_epoch_.*_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')])));
        Nepochs{Count,1} = SubPath{subj}; %#ok<SAGROW>
        Nepochs{Count,2} = regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''); %#ok<SAGROW>
        Nepochs{Count,3} = length(SpkFiles); %#ok<SAGROW>
    end
end
% fprintf('\nEach column corresponds to a different spike type:\n')
% disp([SubPath,Nepochs])
disp(Nepochs)

%% ========== Frequency analysis ==========
FrequencyBands = 1:45;
Delta = 1:4; Theta = 4:8; Alpha = 8:12; Beta = 13:30; Gamma = 30:45;%30:40; initially but bandpass was less narrow in my case
Broad = 1:45;

cowsay('Frequency analysis per epoch...')
for subj = SubjInd
    fprintf('Doing subject %d/%d...\n',subj,length(SubjInd))
    
    % Get .mat file with bad segments (better than .sef file because in
    % .mat bad segments are NaN, whereas they are zeroed in .sef):
    PreprocFile = my_listfiles(fullfile(PreprocPath,SubPath{subj},DataType),'with_bad_TF_i3DS2_SIWMF_avg_ref.mat');
    PreprocFile = char(fullfile(PreprocPath,SubPath{subj},DataType,PreprocFile));
    load(PreprocFile,'ChannelNames')
    
    for spk = 1:length(SpikeLabel{subj})
        
        SpkFiles = cellstr(spm_select('FPList',[fileparts(PreprocFile),filesep,'Epochs'],spm_file(spm_file(PreprocFile,'basename'),'ext','sef','suffix',['_epoch_.*_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')])));
        RndFiles = cellstr(spm_select('FPList',[fileparts(PreprocFile),filesep,'Epochs'],spm_file(spm_file(PreprocFile,'basename'),'ext','sef','suffix',['_random_epoch_.*_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')])));
        
        SpikesEpochs = nan(length(ChannelNames),EpochDuration,length(SpkFiles));
        for ep = 1:length(SpkFiles)
            SpikesEpochs(:,:,ep) = dual_load_sef(SpkFiles{ep});
        end
        RandomsEpochs = nan(length(ChannelNames),EpochDuration,length(SpkFiles));
        for ep = 1:length(RndFiles)
            RandomsEpochs(:,:,ep) = dual_load_sef(RndFiles{ep});
        end
        
        Pxx_pS = nan(length(ChannelNames),length(FrequencyBands),length(SpkFiles));
        Pxx_pR = nan(length(ChannelNames),length(FrequencyBands),length(SpkFiles));
        % Perform frequency analysis with pwelch:
        for ep = 1:size(SpikesEpochs,3)
            prc_for_loop(ep,size(SpikesEpochs,3),1);
            for c = 1:size(SpikesEpochs,1)
%                 fprintf('Doing channel #%d...\n',c)
                [Pxx_p,F_p] = pwelch(vect(SpikesEpochs(c,:,ep)),EpochDuration,0,FrequencyBands,NewFreq);
                Pxx_pS(c,:,ep) = log(Pxx_p); % NB, #RM@FBMlab: TAKING LOG OF FREQUENCY POWER, BECAUSE AVERAGING THE LOG OR CALCULATING THE LOG OF THE AVERAGE IS NOT THE SAME!
                Pxx_p = pwelch(vect(RandomsEpochs(c,:,ep)),EpochDuration,0,FrequencyBands,NewFreq);
                Pxx_pR(c,:,ep) = log(Pxx_p); % NB, #RM@FBMlab: TAKING LOG OF FREQUENCY POWER, BECAUSE AVERAGING THE LOG OR CALCULATING THE LOG OF THE AVERAGE IS NOT THE SAME!
            end
        end
        
        % frequency bands
        Pxx_pS_delta = squeeze(mean(Pxx_pS(:,Delta,:),2));
        Pxx_pR_delta = squeeze(mean(Pxx_pR(:,Delta,:),2));
        Pxx_pS_theta = squeeze(mean(Pxx_pS(:,Theta,:),2));
        Pxx_pR_theta = squeeze(mean(Pxx_pR(:,Theta,:),2));
        Pxx_pS_alpha = squeeze(mean(Pxx_pS(:,Alpha,:),2));
        Pxx_pR_alpha = squeeze(mean(Pxx_pR(:,Alpha,:),2));
        Pxx_pS_beta = squeeze(mean(Pxx_pS(:,Beta,:),2));
        Pxx_pR_beta = squeeze(mean(Pxx_pR(:,Beta,:),2));
        Pxx_pS_gamma = squeeze(mean(Pxx_pS(:,Gamma,:),2));
        Pxx_pR_gamma = squeeze(mean(Pxx_pR(:,Gamma,:),2));
        
        % scaled by broad band power
        Pxx_pS_delta_scaled = Pxx_pS_delta./squeeze(mean(Pxx_pS,2));
        Pxx_pR_delta_scaled = Pxx_pR_delta./squeeze(mean(Pxx_pR,2));
        Pxx_pS_theta_scaled = Pxx_pS_theta./squeeze(mean(Pxx_pS,2));
        Pxx_pR_theta_scaled = Pxx_pR_theta./squeeze(mean(Pxx_pR,2));
        Pxx_pS_alpha_scaled = Pxx_pS_alpha./squeeze(mean(Pxx_pS,2));
        Pxx_pR_alpha_scaled = Pxx_pR_alpha./squeeze(mean(Pxx_pR,2));
        Pxx_pS_beta_scaled = Pxx_pS_beta./squeeze(mean(Pxx_pS,2));
        Pxx_pR_beta_scaled = Pxx_pR_beta./squeeze(mean(Pxx_pR,2));
        Pxx_pS_gamma_scaled = Pxx_pS_gamma./squeeze(mean(Pxx_pS,2));
        Pxx_pR_gamma_scaled = Pxx_pR_gamma./squeeze(mean(Pxx_pR,2));
         
        save(spm_file(PreprocFile,'suffix',['_epochs_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''),'_freq']),'Pxx_pS_*','F_p')
            
        save(spm_file(PreprocFile,'suffix',['_random_epochs_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''),'_freq']),'Pxx_pR_*','F_p')
        
    end
end

%% ========== SEEG frequency analysis ==========
cowsay('SEEG frequency analysis per epoch...')
for subj = SubjInd
    fprintf('Doing subject %d/%d...\n',subj,length(SubjInd))
    
    % Get .mat file with bad segments (better than .sef file because in
    % .mat bad segments are NaN, whereas they are zeroed in .sef):
    SEEG = spm_select('FPList',fullfile(PreprocPath,SubPath{subj},DataType),'^icEEGi_.*down.*\.mat$');
    load(SEEG,'icChannelNames')
    
    % #TODO!!! #TOFIX!!!
    [~,~,icChannelsBip,~,SEEGbip] = bipolar_montage(icChannelNames,2,icEEGd);
    
    for spk = 1:length(SpikeLabel{subj})
        
        SpkFiles = cellstr(spm_select('FPList',[fileparts(SEEG),filesep,'Epochs'],spm_file(spm_file(SEEG,'basename'),'ext','sef','suffix',['_epoch_.*_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')])));
        RndFiles = cellstr(spm_select('FPList',[fileparts(SEEG),filesep,'Epochs'],spm_file(spm_file(SEEG,'basename'),'ext','sef','suffix',['_random_epoch_.*_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|','')])));
        
        SpikesEpochs = nan(length(ChannelNames),EpochDuration,length(SpkFiles));
        for ep = 1:length(SpkFiles)
            SpikesEpochs(:,:,ep) = dual_load_sef(SpkFiles{ep});
        end
        RandomsEpochs = nan(length(ChannelNames),EpochDuration,length(SpkFiles));
        for ep = 1:length(RndFiles)
            RandomsEpochs(:,:,ep) = dual_load_sef(RndFiles{ep});
        end
        
        Pxx_pS = nan(length(ChannelNames),length(FrequencyBands),length(SpkFiles));
        Pxx_pR = nan(length(ChannelNames),length(FrequencyBands),length(SpkFiles));
        % Perform frequency analysis with pwelch:
        for ep = 1:size(SpikesEpochs,3)
            prc_for_loop(ep,size(SpikesEpochs,3),1);
            for c = 1:size(SpikesEpochs,1)
%                 fprintf('Doing channel #%d...\n',c)
                [Pxx_p,F_p] = pwelch(vect(SpikesEpochs(c,:,ep)),EpochDuration,0,FrequencyBands,NewFreq);
                Pxx_pS(c,:,ep) = log(Pxx_p); % NB, #RM@FBMlab: TAKING LOG OF FREQUENCY POWER, BECAUSE AVERAGING THE LOG OR CALCULATING THE LOG OF THE AVERAGE IS NOT THE SAME!
                Pxx_p = pwelch(vect(RandomsEpochs(c,:,ep)),EpochDuration,0,FrequencyBands,NewFreq);
                Pxx_pR(c,:,ep) = log(Pxx_p); % NB, #RM@FBMlab: TAKING LOG OF FREQUENCY POWER, BECAUSE AVERAGING THE LOG OR CALCULATING THE LOG OF THE AVERAGE IS NOT THE SAME!
            end
        end
        
        % frequency bands
        Pxx_pS_delta = squeeze(mean(Pxx_pS(:,Delta,:),2));
        Pxx_pR_delta = squeeze(mean(Pxx_pR(:,Delta,:),2));
        Pxx_pS_theta = squeeze(mean(Pxx_pS(:,Theta,:),2));
        Pxx_pR_theta = squeeze(mean(Pxx_pR(:,Theta,:),2));
        Pxx_pS_alpha = squeeze(mean(Pxx_pS(:,Alpha,:),2));
        Pxx_pR_alpha = squeeze(mean(Pxx_pR(:,Alpha,:),2));
        Pxx_pS_beta = squeeze(mean(Pxx_pS(:,Beta,:),2));
        Pxx_pR_beta = squeeze(mean(Pxx_pR(:,Beta,:),2));
        Pxx_pS_gamma = squeeze(mean(Pxx_pS(:,Gamma,:),2));
        Pxx_pR_gamma = squeeze(mean(Pxx_pR(:,Gamma,:),2));
        
        % scaled by broad band power
        Pxx_pS_delta_scaled = Pxx_pS_delta./squeeze(mean(Pxx_pS,2));
        Pxx_pR_delta_scaled = Pxx_pR_delta./squeeze(mean(Pxx_pR,2));
        Pxx_pS_theta_scaled = Pxx_pS_theta./squeeze(mean(Pxx_pS,2));
        Pxx_pR_theta_scaled = Pxx_pR_theta./squeeze(mean(Pxx_pR,2));
        Pxx_pS_alpha_scaled = Pxx_pS_alpha./squeeze(mean(Pxx_pS,2));
        Pxx_pR_alpha_scaled = Pxx_pR_alpha./squeeze(mean(Pxx_pR,2));
        Pxx_pS_beta_scaled = Pxx_pS_beta./squeeze(mean(Pxx_pS,2));
        Pxx_pR_beta_scaled = Pxx_pR_beta./squeeze(mean(Pxx_pR,2));
        Pxx_pS_gamma_scaled = Pxx_pS_gamma./squeeze(mean(Pxx_pS,2));
        Pxx_pR_gamma_scaled = Pxx_pR_gamma./squeeze(mean(Pxx_pR,2));
         
        save(spm_file(SEEG,'suffix',['_epochs_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''),'_freq']),'Pxx_pS_*','F_p')
            
        save(spm_file(SEEG,'suffix',['_random_epochs_for_',regexprep(SpikeLabel{subj}{spk},'\$|\^|"|\|',''),'_freq']),'Pxx_pR_*','F_p')
        
    end
end

%% D. Brunet's pipeline as of Spring 2019
% 1. filtering
% 2. downsampling
% 3. bad epochs (crazy electrodes)
% 4. ICA  (MNE)
% 5. interpolation
% 6. spatial filter
% 7. re-referencing
% 8. Scanning for Bad Epochs (Cartool): eye blinks, potentially other things
