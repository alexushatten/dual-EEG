
SPIfilePaths = {'E:\FS_subjects_DONE\sub-50\LSMAC_no_cereb\More\T1_fixed.spi'
                'E:\FS_subjects_DONE\sub-48\LSMAC_no_cereb\More\T1_fixed.spi'
                'E:\FS_subjects_DONE\sub-82\LSMAC_no_cereb\More\T1_fixed.spi'
                'E:\FS_subjects_DONE\sub-01\LSMAC_no_cereb\More\T1_no_headphones_masked.spi'
                'E:\FS_subjects_DONE\sub-35\LSMAC_no_cereb\More\T1_no_headphones.spi'
                'E:\FS_subjects_DONE\sub-96\LSMAC_no_cereb\More\T1_fixed.spi'
                'E:\FS_subjects_DONE\sub-33\LSMAC_no_cereb\More\T1.spi'
                'E:\FS_subjects_DONE\sub-11\LSMAC_no_cereb\More\T1.spi'
                'E:\FS_subjects_DONE\sub-34\LSMAC_no_cereb\More\T1.spi'
                'E:\FS_subjects_DONE\sub-12\LSMAC_no_cereb\More\T1.spi'};

ISfilePaths = {'E:\FS_subjects_DONE\sub-50\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-48\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-82\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-01\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-35\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-96\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-33\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-11\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-34\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'
               'E:\FS_subjects_DONE\sub-12\LSMAC_no_cereb\LSMAC_no_cereb.Laura.is'};

ELSfilePaths = {'E:\FS_subjects_DONE\sub-50\elec_recon\sub50_align_equ_de174.els'
                'E:\FS_subjects_DONE\sub-48\elec_recon\sub48_align_equ_de222.els'
                'E:\FS_subjects_DONE\sub-82\elec_recon\sub82_align_equ_de160_fix.els'
                'E:\FS_subjects_DONE\sub-01\elec_recon\sub01_align_equ_de81.els'
                'E:\FS_subjects_DONE\sub-35\elec_recon\sub35_align_equ_de112.els'
                'E:\FS_subjects_DONE\sub-96\elec_recon\sub96_align_equ_de194.els'
                'E:\FS_subjects_DONE\sub-33\elec_recon\sub33_align_equ_de128.els'
                'E:\FS_subjects_DONE\sub-11\elec_recon\sub11_align_equ_de81.els'
                'E:\FS_subjects_DONE\sub-34\elec_recon\sub34_align_equ_de101.els'
                'E:\FS_subjects_DONE\sub-12\elec_recon\sub12_align_equ_de49.els'};

AltSPIfilePaths = {'E:\FS_subjects_DONE\sub-50\LSMAC_cereb\More\T1_fixed.spi'
                   'E:\FS_subjects_DONE\sub-48\LSMAC_cereb\More\T1_fixed.spi'
                   'E:\FS_subjects_DONE\sub-82\LSMAC_cereb\More\T1_fixed.spi'
                   'E:\FS_subjects_DONE\sub-01\LSMAC_cereb\More\T1_no_headphones_masked.spi'
                   'E:\FS_subjects_DONE\sub-35\LSMAC_cereb\More\T1_no_headphones.spi'
                   'E:\FS_subjects_DONE\sub-96\LSMAC_cereb\More\T1_fixed.spi'
                   'E:\FS_subjects_DONE\sub-33\LSMAC_cereb\More\T1.spi'
                   'E:\FS_subjects_DONE\sub-11\LSMAC_cereb\More\T1.spi'
                   'E:\FS_subjects_DONE\sub-34\LSMAC_cereb\More\T1.spi'
                   'E:\FS_subjects_DONE\sub-12\LSMAC_cereb\More\T1.spi'};

AltISfilePaths = {'E:\FS_subjects_DONE\sub-50\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-48\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-82\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-01\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-35\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-96\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-33\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-11\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-34\LSMAC_cereb\LSMAC_cereb.Laura.is'
                  'E:\FS_subjects_DONE\sub-12\LSMAC_cereb\LSMAC_cereb.Laura.is'};

regexprep(ISfilePaths,'Laura.is','Loreta.is') % To use LORETA instead of LAURA inverse...

%% Definitions

Pre = 199;
Post = 200;

%% File paths
ProcessedHDeegFilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^hdEEG_.*down200hz.*i3DS2_SIWMF_avg_ref.*\.mat',1,...
    'basename','(freq|ICA_infos|artefacts|epochs|hilbert)','out'));
ArtefactsInfoPaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^hdEEG_.*resting_artefacts.*\.mat',1,...
    'basename','(freq|ICA_infos|epochs|hilbert)','out'));
SEEGfilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^icEEGi_.*down200hz.*.*\.mat',1,...
    'basename','(freq|ICA_infos|artefacts|epochs|hilbert)','out'));
EpochingFilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^hdEEG_.*down200hz.*i3DS2_SIWMF.*epochs.*\.mat',1,...
    'basename','(freq|ICA_infos|artefacts|hilbert)','out'));
PatFol = cellfun(@spm_file,cellfun(@fileparts,cellfun(@fileparts,cellfun(@fileparts,EpochingFilePaths,'uniformoutput',0),'uniformoutput',0),'uniformoutput',0),repmat({'basename'},numel(EpochingFilePaths),1),'uniformoutput',0);
GroupByPat = match_vectors(unique(PatFol),PatFol,0);
EpochingFilePaths_BKP = EpochingFilePaths; EpochingFilePaths = {};
for pat = 1:length(GroupByPat)
    EpochingFilePaths{pat} = EpochingFilePaths_BKP(GroupByPat{pat}); %#ok<SAGROW>
end; EpochingFilePaths = EpochingFilePaths';

%% Processing loop

% % Instead of:
% [Pxx_p,F_p] = pwelch(vect(Sources(1,:,1)),200,0,1:45,200);
% % Use FieldTrip for all frequency and time-frequency analyses

% NB: the number of channels of SEEG files might be different from the
% number of channels in the .els file, because MKR channels were kept in
% SEEG files => discard them if using get_max_chan_TF.m !

% NB: we will use convolutions with Morlet wavelets, because Hilbert
% transform on bandpassed signals is not recommended for short recordings
% such as epochs (https://neuroimage.usc.edu/brainstorm/Tutorials/TimeFrequency).

% =====================================================================
% ============== ESI using 1st eigenvariate (SVD method) ==============
% ==============    but based on L-corner of Tikhonov    ==============
% ==============    regularization values against the    ==============
% ==============   average norm of all solution points   ==============
% =====================================================================

% ==================================================================== %
% ========== Here we use the head model without cerebellum =========== %
% ==================================================================== %

for pat = 1:length(ProcessedHDeegFilePaths)
    fprintf('Doing patient %s...\n',num2str(pat))
    load(ProcessedHDeegFilePaths{pat})
    load(ArtefactsInfoPaths{pat})
    clear Epochs
    for spk = 1:length(EpochingFilePaths{pat})
        fprintf('Doing condition %s...\n',num2str(spk))
        load(EpochingFilePaths{pat}{spk})
        SamplingRate = NewFreq;
        
        FileBasename = spm_file(EpochingFilePaths{pat}{spk},'basename');
        
        if strcmpi(FileBasename(end-12:end),'_spike_epochs')
            Epochs = get_EEG_epochs(EEGfull,SpkStart,Pre,Post);
            EpochStart = SpkStart-Pre;
            EpochEnd = SpkStart+Post;
            EpochTiming = SpikesT;
        elseif strcmpi(FileBasename(end-12:end),'random_epochs')
            Epochs = get_EEG_epochs(EEGfull,RandomStart,Pre,Post);
            EpochStart = RandomStart-Pre;
            EpochEnd = RandomStart+Post;
            EpochTiming = RandomsT;
        end
        
        % Pre-allocate Sources:
        IS = fbmlab_read_is(ISfilePaths{pat});
        Sources = nan(IS.numsolutionpoints,size(Epochs,2),size(Epochs,3));
        % Estimates sources for all epochs:
        fprintf('Doing sources for epochs...\n')
        parfor ep = 1:size(Epochs,3)
%         for ep = 1:size(Epochs,3)
%             prc_for_loop(ep,size(Epochs,3),1);
            [~,~,RISsvd,OptimalReg] = compute_inverse(Epochs(:,:,ep),ISfilePaths{pat},'optimal');
            Sources(:,:,ep) = RISsvd{OptimalReg};
        end
        
        if strcmpi(FileBasename(end-12:end),'_spike_epochs')
            OutputFilename = regexprep(regexprep(regexprep(EpochingFilePaths{pat}{spk},'_spike_epochs.mat','.mat'),'hdEEG_','ESI_no_cereb_'),'_with_bad_TF','');
        elseif strcmpi(FileBasename(end-12:end),'random_epochs')
            OutputFilename = regexprep(regexprep(regexprep(EpochingFilePaths{pat}{spk},'_random_epochs.mat','.mat'),'hdEEG_','ESI_no_cereb_'),'_with_bad_TF','');
        end
        
        fprintf('Saving output...\n')
        save(OutputFilename,'Sources','-v7.3')
        
    end
end

% ... Now, the same but with ROI space (svd again but from SPs to ROIs)...

% ===============================================================
% ===========       Reconstruct signal in ROIs        ===========
% ==== (running SVD on vectorized orientations & regions is  ====
% ====  the same as running SVD on regions where SVD already ====
% ====  extracted the 1st eigenvector across orientations)   ====
% ===============================================================

for f = 1:length(ESIfilePaths)
    fprintf('Doing patient %s...\n',num2str(pat))
    
    % ====== Load sources & ROIs ======
    load(ESIfilePaths{f})
    
    % Pre-allocate:
    ROIs = read_rois_Cartool(spm_file(SPIfilePaths{WhichSubj(f)},'ext','rois'));
    SourcesROIs = nan(ROIs.NumberOfROIs,size(Sources,2),size(Sources,3));
    
    fprintf('Computing sources in ROIs...\n')
    parfor ep = 1:size(Sources,3)
        SourcesROIs(:,:,ep) = ESI2ROIs(Sources(:,:,ep),spm_file(SPIfilePaths{WhichSubj(f)},'ext','rois'));
    end
    
    fprintf('Saving output...\n')
    save(regexprep(ESIfilePaths{f},'ESI_no_cereb_','ESI_no_cereb_ROIs_'),'SourcesROIs','-v7.3')
    
end

ROIsESIfilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^ESI_.*ROIs_.*down200hz.*epochs.*\.mat',1,...
    'basename','(freq|ICA_infos|artefacts)','out'));
ROIsNames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71';'72';'73';'74';'75';'76';'77';'78';'79';'80';'81';'82'};

%======================================================================%
% Let's redo that on continuous EEG to perform time-frequency analyses %
%======================================================================%

for pat = 1:length(ProcessedHDeegFilePaths)
    fprintf('Doing patient %s...\n',num2str(pat))
    load(ProcessedHDeegFilePaths{pat})
    EEGnoNaN = EEGfull(:,any(~isnan(EEGfull)));
    [ContinuousSourcesROIs,OptimalReg] = ...
        compute_inverse_in_ROIs(EEGnoNaN,ISfilePaths{pat},spm_file(SPIfilePaths{pat},'ext','rois'));
    OutputFilename = regexprep(ProcessedHDeegFilePaths{pat},'hdEEG_','ESI_no_cereb_ROIs_');
    save(OutputFilename,'ContinuousSourcesROIs','OptimalReg','-v7.3')
end

ContROIsESIfilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^ESI_.*ROIs_.*down200hz.*_avg_ref\.mat',1,...
    'basename','(freq|ICA_infos|artefacts)','out'),...
    'path','Epochs','out');

% ==================================================================== %
% .......... Now let's use FieldTrip instead of pwelch to ............ %
% ......... perform all subsequent frequency and TF analyses ......... %
% ==================================================================== %
% ============== Frequency and time-frequency analyses =============== %
% ==================================================================== %

addpath('E:\code\MATLAB\MATLAB_toolboxes\fieldtrip_20170612'); ft_defaults;

ProcessedICeegFilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^icEEGi_.*down200hz.*\.mat',1,...
    'basename','(freq|ICA_infos|artefacts|epochs|hilbert)','out'));
EpochsICeegFilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^icEEGi_.*down200hz.*epochs.*\.mat',1,...
    'basename','(freq|ICA_infos|artefacts)','out'));
ESIfilePaths = cellstr(my_super_listfiles('E:\dualEEG\3_Preprocessed',...
    '^ESI_.*down200hz.*epochs.*\.mat',1,...
    'basename','(freq|ICA_infos|artefacts|ROI)','out'));

SubjPathParsed = get_paths_part(EpochsICeegFilePaths,4);
WhichSubj = match_vectors(SubjPathParsed,unique(SubjPathParsed),1); % this works because SubjPathParsed is already sorted

FrequencyBands = 1:45;
Delta = 1:4; Theta = 4:8; Alpha = 8:12; Beta = 13:30; Gamma = 30:45;%31:40;%30:45;%30:40; initially but bandpass was less narrow in my case
% Broad = 1:45;
% Edits based on http://www.scholarpedia.org/article/Electroencephalogram
% Another resource: https://www.medicine.mcgill.ca/physio/vlab/biomed_signals/eeg_n.htm
% Overlap between frequency bands was initially removed, but there were
% overlaps across multiple bands from the definition anyway...

% It is annoying to do it everytime per subject and per condition, we will
% loop over files. This is ok since the order of other paths is
% alphabetical!

for f = 1:length(ESIfilePaths)
    
    % ====== Load sources ======
    load(ESIfilePaths{f})
    % ====== Load SEEG ======
    load(EpochsICeegFilePaths{f})
    % ====== Load SEEG channels info ======
    load(ProcessedICeegFilePaths{WhichSubj(f)},'icChannelNames')
    % ====== Make bipolar montage (without bad channels) ======
    SEEGchannels = icChannelNames(cellfun(@isempty,regexp(icChannelNames,'Mkr')));
    SEEGchannels = SEEGchannels(cellfun(@isempty,regexp(SEEGchannels,'ecg')));
    
    % Pre-allocate:
    Temp = bipolar_montage(SEEGchannels,2);
    bicEEGepochs = nan(numel(Temp),size(icEEGepochs,2),size(icEEGepochs,3));
    for ep = 1:size(icEEGepochs,3)
        [~,~,bip_labels,bip_labels2,bicEEGepochs(:,:,ep)] = bipolar_montage(SEEGchannels,2,icEEGepochs(:,:,ep));
    end
    
    % ===============================================================
    % ===========        Frequency analysis on ESI        ===========
    % ===============================================================
    
    % =========== Re-format data to match FieldTrip input ===========
    Hdr.Fs = 200;
    Hdr.nChans = size(Sources,1);
    Hdr.nSamplesPre = 0;
    Hdr.nSamples = numel(Sources);
    Hdr.nTrials = size(Sources,3);
    Hdr.label = regexprep(cellstr(spm_file(num2str(vect(1:size(Sources,1))),'prefix','SP_')),' ','');
    Hdr.orig = []; % won't be used by ft_freqanalysis anyway...
    Hdr.chantype = cellstr(repmat('source',size(Sources,1),1)); % will be used by ft_freqanalysis indirectly, via call to ft_datatype and getdimord !
    Hdr.chanunit = cellstr(repmat('uV',size(Sources,1),1)); % won't be used by ft_freqanalysis anyway...
    
    Data = [];
    Data.hdr = Hdr;
    Data.label = Hdr.label;
    % Data.time = (0:(Data.hdr.nSamples-1))/Data.hdr.Fs; % => NO, wrong!
    % Data.time = (0:(size(Sources,2)-1))/Data.hdr.Fs; % => NO, should match structure of "trial" field!
    Count = 0;
    for t = 1:Hdr.nTrials
        Data.time{t} = (Count:(Count+size(Sources,2)-1))/Data.hdr.Fs;
        Count = size(Data.time{t});
    end
    for t = 1:Hdr.nTrials
        Data.trial{t} = Sources(:,:,t);
    end
    % Data.trial = Sources;
    Data.fsample = Hdr.Fs;
    Data.sampleinfo = [vect(1:size(Sources,2):size(Sources,2)*size(Sources,3)),vect(size(Sources,2):size(Sources,2):size(Sources,2)*size(Sources,3))];
    % Data.trialinfo = [vect(1:size(Sources,2):size(Sources,2)*size(Sources,3)),vect(size(Sources,2):size(Sources,2):size(Sources,2)*size(Sources,3))];
    
    % The datatype has to be defined, otherwise ft_freqanalysis.m won't let us
    % perform the task...
    SPxyz = read_spi_Cartool(SPIfilePaths{WhichSubj(f)});
    
    Data.pos = [SPxyz.x',SPxyz.y',SPxyz.z'];
    % Data.ori = 1; % FieldTrip assumes number of orientations (of dipoles) to
    % % be 3 by default, but does it solve the dimord issue with trial?... Try again!
    
    % Data.dimord = 'rpt_chan_time'; % trials (repetitions) x channels x time
    % Data.dimord = 'rpt_pos_time'; % trials (repetitions) x channels x time
    Data.trialdimord = 'rpt_pos_time';
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.foi = FrequencyBands;
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 2;
    cfg.channel = Data.label;
    cfg.keeptrials = 'yes';
    ESIfreqResults = ft_freqanalysis(cfg,Data);
    
    % Average over frequency bands corresponding to delta,
    % theta, alpha, beta and gamma, but log-transform before:
    ESIdeltaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Delta)),3);
    ESIthetaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Theta)),3);
    ESIalphaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Alpha)),3);
    ESIbetaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Beta)),3);
    ESIgammaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Gamma)),3);
    
%     % ===============================================================
%     % ===========      Time-frequency analysis on ESI     ===========
%     % ===============================================================
%     cfg = [];
%     cfg.method = 'wavelet';
%     cfg.output = 'pow';
%     cfg.foi = FrequencyBands;
%     cfg.width = 4;
%     cfg.toi = 'all';
%     cfg.channel = Data.label;
%     cfg.keeptrials = 'yes';
%     ESItimeFreqResults = ft_freqanalysis(cfg,Data);
%     
%     % Average over frequency bands corresponding to delta,
%     % theta, alpha, beta and gamma, but log-transform before:
%     ESIdeltaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Delta,:)),3));
%     ESIthetaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Theta,:)),3));
%     ESIalphaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Alpha,:)),3));
%     ESIbetaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Beta,:)),3));
%     ESIgammaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Gamma,:)),3));
    
    % ===============================================================
    % ===========       Frequency analysis on SEEG        ===========
    % ===============================================================
    Hdr.Fs = 200;
    Hdr.nChans = size(bicEEGepochs,1);
    Hdr.nSamplesPre = 0;
    Hdr.nSamples = numel(bicEEGepochs);
    Hdr.nTrials = size(bicEEGepochs,3);
    Hdr.label = bip_labels';
    Hdr.orig = []; % won't be used by ft_freqanalysis anyway...
    Hdr.chantype = cellstr(repmat('seeg',size(bicEEGepochs,1),1)); % will be used by ft_freqanalysis indirectly, via call to ft_datatype and getdimord !
    Hdr.chanunit = cellstr(repmat('uV',size(bicEEGepochs,1),1)); % won't be used by ft_freqanalysis anyway...
    
    Data = [];
    Data.hdr = Hdr;
    Data.label = Hdr.label;
    % Data.time = (0:(Data.hdr.nSamples-1))/Data.hdr.Fs; % => NO, wrong!
    % Data.time = (0:(size(Sources,2)-1))/Data.hdr.Fs; % => NO, should match structure of "trial" field!
    Count = 0;
    for t = 1:Hdr.nTrials
        Data.time{t} = (Count:(Count+size(bicEEGepochs,2)-1))/Data.hdr.Fs;
        Count = size(Data.time{t});
    end
    for t = 1:Hdr.nTrials
        Data.trial{t} = bicEEGepochs(:,:,t);
    end
    % Data.trial = Sources;
    Data.fsample = Hdr.Fs;
    Data.sampleinfo = [vect(1:size(bicEEGepochs,2):size(bicEEGepochs,2)*size(bicEEGepochs,3)),vect(size(bicEEGepochs,2):size(bicEEGepochs,2):size(bicEEGepochs,2)*size(bicEEGepochs,3))];
    
    Data.trialdimord = 'rpt_pos_time';
    
    [x,y,z,name,ClusterName,Nclus,FullName] = read_els_file(ELSfilePaths{WhichSubj(f)});
    
    ChanIdx = reshape(match_vectors(bip_labels2(:),FullName(:),0),size(bip_labels2));
    
    LHSchan = get_chan_coordinates(ELSfilePaths{WhichSubj(f)},ChanIdx(:,1));
    RHSchan = get_chan_coordinates(ELSfilePaths{WhichSubj(f)},ChanIdx(:,2));
    
    SEEGcoor = wmean_mat(LHSchan,RHSchan,[1 1]);
    Data.pos = SEEGcoor;
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.foi = FrequencyBands;
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 2;
    cfg.channel = Data.label;
    cfg.keeptrials = 'yes';
    SEEGfreqResults = ft_freqanalysis(cfg,Data);
    
    % Average over frequency bands corresponding to delta,
    % theta, alpha, beta and gamma, but log-transform before:
    SEEGdeltaFreqResults = nanmean(log(SEEGfreqResults.powspctrm(:,:,Delta)),3);
    SEEGthetaFreqResults = nanmean(log(SEEGfreqResults.powspctrm(:,:,Theta)),3);
    SEEGalphaFreqResults = nanmean(log(SEEGfreqResults.powspctrm(:,:,Alpha)),3);
    SEEGbetaFreqResults = nanmean(log(SEEGfreqResults.powspctrm(:,:,Beta)),3);
    SEEGgammaFreqResults = nanmean(log(SEEGfreqResults.powspctrm(:,:,Gamma)),3);
    
%     % ===============================================================
%     % ===========     Time-frequency analysis on SEEG     ===========
%     % ===============================================================
%     cfg = [];
%     cfg.method = 'wavelet';
%     cfg.output = 'pow';
%     cfg.foi = FrequencyBands;
%     cfg.width = 4;
%     cfg.toi = 'all';
%     cfg.channel = Data.label;
%     cfg.keeptrials = 'yes';
%     SEEGtimeFreqResults = ft_freqanalysis(cfg,Data);
%     
%     % Average over frequency bands corresponding to delta,
%     % theta, alpha, beta and gamma, but log-transform before:
%     SEEGdeltaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,:,Delta,:)),3));
%     SEEGthetaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,:,Theta,:)),3));
%     SEEGalphaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,:,Alpha,:)),3));
%     SEEGbetaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,:,Beta,:)),3));
%     SEEGgammaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,:,Gamma,:)),3));
    
    % Save output:
    OutputFilename = regexprep(regexprep(ESIfilePaths{f},'ESI_',''),'_ICA_recon_i3DS2_SIWMF_avg_ref','_freq_SPs_SEEG');
    
    save(OutputFilename,...
    'ESIfreqResults','ESIdeltaFreqResults','ESIthetaFreqResults','ESIalphaFreqResults','ESIbetaFreqResults','ESIgammaFreqResults',...
    'SEEGfreqResults','SEEGdeltaFreqResults','SEEGthetaFreqResults','SEEGalphaFreqResults','SEEGbetaFreqResults','SEEGgammaFreqResults',...
    '-v7.3');
%     'ESItimeFreqResults','ESIdeltaTimeFreqResults','ESIthetaTimeFreqResults','ESIalphaTimeFreqResults','ESIbetaTimeFreqResults','ESIgammaTimeFreqResults',...
%     'SEEGtimeFreqResults','SEEGdeltaTimeFreqResults','SEEGthetaTimeFreqResults','SEEGalphaTimeFreqResults','SEEGbetaTimeFreqResults','SEEGgammaTimeFreqResults',...
    
end


for f = 1:length(ROIsESIfilePaths)
    
    % ====== Load sources in ROIs ======
    load(ROIsESIfilePaths{f})

    % ===============================================================
    % ===========        Frequency analysis on ESI        ===========
    % ===============================================================
    
    % =========== Re-format data to match FieldTrip input ===========
    Hdr.Fs = 200;
    Hdr.nChans = size(SourcesROIs,1);
    Hdr.nSamplesPre = 0;
    Hdr.nSamples = numel(SourcesROIs);
    Hdr.nTrials = size(SourcesROIs,3);
    Hdr.label = ROIsNames;
    Hdr.orig = []; % won't be used by ft_freqanalysis anyway...
    Hdr.chantype = cellstr(repmat('source',size(SourcesROIs,1),1)); % will be used by ft_freqanalysis indirectly, via call to ft_datatype and getdimord !
    Hdr.chanunit = cellstr(repmat('uV',size(SourcesROIs,1),1)); % won't be used by ft_freqanalysis anyway...
    
    Data = [];
    Data.hdr = Hdr;
    Data.label = Hdr.label;
    % Data.time = (0:(Data.hdr.nSamples-1))/Data.hdr.Fs; % => NO, wrong!
    % Data.time = (0:(size(SourcesROIs,2)-1))/Data.hdr.Fs; % => NO, should match structure of "trial" field!
    Count = 0;
    for t = 1:Hdr.nTrials
        Data.time{t} = (Count:(Count+size(SourcesROIs,2)-1))/Data.hdr.Fs;
        Count = size(Data.time{t});
    end
    for t = 1:Hdr.nTrials
        Data.trial{t} = SourcesROIs(:,:,t);
    end
    % Data.trial = SourcesROIs;
    Data.fsample = Hdr.Fs;
    Data.sampleinfo = [vect(1:size(SourcesROIs,2):size(SourcesROIs,2)*size(SourcesROIs,3)),vect(size(SourcesROIs,2):size(SourcesROIs,2):size(SourcesROIs,2)*size(SourcesROIs,3))];
    % Data.trialinfo = [vect(1:size(SourcesROIs,2):size(SourcesROIs,2)*size(SourcesROIs,3)),vect(size(SourcesROIs,2):size(SourcesROIs,2):size(SourcesROIs,2)*size(SourcesROIs,3))];
    
    % The datatype has to be defined, otherwise ft_freqanalysis.m won't let us
    % perform the task...
    SPxyz = read_spi_Cartool(regexprep(SPIfilePaths{WhichSubj(f)},'.spi','.rois.spi'));
    
    Data.pos = [SPxyz.x',SPxyz.y',SPxyz.z'];
    % Data.ori = 1; % FieldTrip assumes number of orientations (of dipoles) to
    % % be 3 by default, but does it solve the dimord issue with trial?... Try again!
    
    % Data.dimord = 'rpt_chan_time'; % trials (repetitions) x channels x time
    % Data.dimord = 'rpt_pos_time'; % trials (repetitions) x channels x time
    Data.trialdimord = 'rpt_pos_time';
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.foi = FrequencyBands;
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 2;
    cfg.channel = Data.label;
    cfg.keeptrials = 'yes';
    ESIfreqResults = ft_freqanalysis(cfg,Data);
    
    % Average over frequency bands corresponding to delta,
    % theta, alpha, beta and gamma, but log-transform before:
    ESIdeltaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Delta)),3);
    ESIthetaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Theta)),3);
    ESIalphaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Alpha)),3);
    ESIbetaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Beta)),3);
    ESIgammaFreqResults = nanmean(log(ESIfreqResults.powspctrm(:,:,Gamma)),3);
    
%     % ===============================================================
%     % ===========      Time-frequency analysis on ESI     ===========
%     % ===============================================================
%     cfg = [];
%     cfg.method = 'wavelet';
%     cfg.output = 'pow';
%     cfg.foi = FrequencyBands;
%     cfg.width = 4;
%     cfg.toi = 'all';
%     cfg.channel = Data.label;
%     cfg.keeptrials = 'yes';
%     ESItimeFreqResults = ft_freqanalysis(cfg,Data);
%     
%     % Average over frequency bands corresponding to delta,
%     % theta, alpha, beta and gamma, but log-transform before:
%     ESIdeltaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Delta,:)),3));
%     ESIthetaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Theta,:)),3));
%     ESIalphaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Alpha,:)),3));
%     ESIbetaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Beta,:)),3));
%     ESIgammaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,:,Gamma,:)),3));
    
    % Save output:
    OutputFilename = regexprep(regexprep(ROIsESIfilePaths{f},'ESI_ROIs_',''),'_ICA_recon_i3DS2_SIWMF_avg_ref','_freq_ROIs');
    
    save(OutputFilename,...
    'ESIfreqResults','ESIdeltaFreqResults','ESIthetaFreqResults','ESIalphaFreqResults','ESIbetaFreqResults','ESIgammaFreqResults',...
    '-v7.3');
%     'ESItimeFreqResults','ESIdeltaTimeFreqResults','ESIthetaTimeFreqResults','ESIalphaTimeFreqResults','ESIbetaTimeFreqResults','ESIgammaTimeFreqResults',...
    
end

% The results above do not estimate properly power in low-frequency bands...
%% Alternative, assuming I have continuous signals (reconstruction of whole time courses needed, will be memory intensive) %%
% ...apply hilbert transform on band-pass filtered signals:
Test_1_5 = ft_preproc_bandpassfilter(bsxfun(@minus,ContinuousSourcesROIs,mean(ContinuousSourcesROIs,2)),200,[1 5],4,'but','twopass');
% Then the magnitude gives the power of the frequency (we've filtered at):
FreqPower_1_5 = abs(hilbert(Test_1_5(1,:)));

for f = 1:length(ContROIsESIfilePaths)
    
    % ====== Load sources ======
    load(ContROIsESIfilePaths{f})
    % ====== Load SEEG ======
    load(ProcessedICeegFilePaths{f})
    % ====== Make bipolar montage (without bad channels) ======
    SEEGchannels = icChannelNames(cellfun(@isempty,regexp(icChannelNames,'Mkr')));
    SEEGchannels = SEEGchannels(cellfun(@isempty,regexp(SEEGchannels,'ecg')));
    icEEGnoNaN = icEEGd(:,any(~isnan(icEEGd)));
    [~,~,bip_labels,bip_labels2,bipSEEG] = bipolar_montage(SEEGchannels,2,icEEGnoNaN);
    
    % ===============================================================
    % ===========      Time-frequency analysis on ESI     ===========
    % ===============================================================
    
    % =========== Re-format data to match FieldTrip input ===========
    Hdr = [];
    Hdr.Fs = 200;
    Hdr.nChans = size(ContinuousSourcesROIs,1);
    Hdr.nSamplesPre = 0;
    Hdr.nSamples = size(ContinuousSourcesROIs,2);
    Hdr.nTrials = 1;
    Hdr.label = regexprep(cellstr(num2str(vect(1:size(ContinuousSourcesROIs,1)))),' ','');
    Hdr.orig = []; % won't be used by ft_freqanalysis anyway...
    Hdr.chantype = cellstr(repmat('source',size(ContinuousSourcesROIs,1),1)); % will be used by ft_freqanalysis indirectly, via call to ft_datatype and getdimord !
    Hdr.chanunit = cellstr(repmat('uV',size(ContinuousSourcesROIs,1),1)); % won't be used by ft_freqanalysis anyway...
    
    Data = [];
    Data.hdr = Hdr;
    Data.label = Hdr.label;
    Data.time{1} = (0:(size(ContinuousSourcesROIs,2)-1))/Data.hdr.Fs;
    Data.trial{1} = ContinuousSourcesROIs;
    Data.fsample = Hdr.Fs;
    Data.sampleinfo = [1 size(ContinuousSourcesROIs,2)];
    
    SPxyz = read_spi_Cartool(SPIfilePaths{f});
    
    Data.pos = [SPxyz.x',SPxyz.y',SPxyz.z'];
    
    % Data.dimord = 'rpt_chan_time'; % trials (repetitions) x channels x time
    % Data.dimord = 'rpt_pos_time'; % trials (repetitions) x channels x time
    Data.trialdimord = 'rpt_pos_time';
    
    cfg = [];
    cfg.method = 'wavelet';
    cfg.output = 'pow';
    cfg.foi = FrequencyBands;
    cfg.width = 4;
    cfg.toi = 'all';
    cfg.channel = Data.label;
%     cfg.keeptrials = 'yes';
    ESItimeFreqResults = ft_freqanalysis(cfg,Data);
    
    % Average over frequency bands corresponding to delta,
    % theta, alpha, beta and gamma, but log-transform before:
    ESIdeltaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,Delta,:)),2));
    ESIthetaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,Theta,:)),2));
    ESIalphaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,Alpha,:)),2));
    ESIbetaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,Beta,:)),2));
    ESIgammaTimeFreqResults = squeeze(nanmean(log(ESItimeFreqResults.powspctrm(:,Gamma,:)),2));
    
    % ===============================================================
    % ===========     Time-frequency analysis on SEEG     ===========
    % ===============================================================

    Hdr.Fs = 200;
    Hdr.nChans = size(bipSEEG,1);
    Hdr.nSamplesPre = 0;
    Hdr.nSamples = size(bipSEEG,2);
    Hdr.nTrials = 1;
    Hdr.label = bip_labels';
    Hdr.orig = []; % won't be used by ft_freqanalysis anyway...
    Hdr.chantype = cellstr(repmat('seeg',size(bipSEEG,1),1)); % will be used by ft_freqanalysis indirectly, via call to ft_datatype and getdimord !
    Hdr.chanunit = cellstr(repmat('uV',size(bipSEEG,1),1)); % won't be used by ft_freqanalysis anyway...
    
    Data = [];
    Data.hdr = Hdr;
    Data.label = Hdr.label;
    Data.time{1} = (0:(size(bipSEEG,2)-1))/Data.hdr.Fs;
    Data.trial{1} = bipSEEG;
    Data.fsample = Hdr.Fs;
    Data.sampleinfo = [1 size(bipSEEG,2)];
    
    Data.trialdimord = 'rpt_pos_time';
    
    [x,y,z,name,ClusterName,Nclus,FullName] = read_els_file(ELSfilePaths{f});
    
    ChanIdx = reshape(match_vectors(bip_labels2(:),FullName(:),0),size(bip_labels2));
    
    LHSchan = get_chan_coordinates(ELSfilePaths{WhichSubj(f)},ChanIdx(:,1));
    RHSchan = get_chan_coordinates(ELSfilePaths{WhichSubj(f)},ChanIdx(:,2));
    
    SEEGcoor = wmean_mat(LHSchan,RHSchan,[1 1]);
    Data.pos = SEEGcoor;
    
    cfg = [];
    cfg.method = 'wavelet';
    cfg.output = 'pow';
    cfg.foi = FrequencyBands;
    cfg.width = 4;
    cfg.toi = 'all';
    cfg.channel = Data.label;
%     cfg.keeptrials = 'yes';
    SEEGtimeFreqResults = ft_freqanalysis(cfg,Data);
    
    % Average over frequency bands corresponding to delta,
    % theta, alpha, beta and gamma, but log-transform before:
    SEEGdeltaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,Delta,:)),2));
    SEEGthetaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,Theta,:)),2));
    SEEGalphaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,Alpha,:)),2));
    SEEGbetaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,Beta,:)),2));
    SEEGgammaTimeFreqResults = squeeze(nanmean(log(SEEGtimeFreqResults.powspctrm(:,Gamma,:)),2));
    
    % Save output:
    OutputFilename = regexprep(regexprep(regexprep(ContROIsESIfilePaths,'_ICA_recon_with_bad_TF_i3DS2_SIWMF_avg_ref','_timefreq'),'ESI_no_cereb_ROIs_','ESI_no_cereb_ROIs_SEEG_'),'_with_bad_TF','');
    
    save(OutputFilename,...
    'ESItimeFreqResults','ESIdeltaTimeFreqResults','ESIthetaTimeFreqResults','ESIalphaTimeFreqResults','ESIbetaTimeFreqResults','ESIgammaTimeFreqResults',...
    'SEEGtimeFreqResults','SEEGdeltaTimeFreqResults','SEEGthetaTimeFreqResults','SEEGalphaTimeFreqResults','SEEGbetaTimeFreqResults','SEEGgammaTimeFreqResults',...
    '-v7.3');
    
end

% ~~~~~~~~~~~ Inverses computed in Cartool are corrupted because of a bug ~~~~~~~~~~~ %
% ~~~~~~~~~~~~~~~~~~~~~ Let's do it in sensor space, instead... ~~~~~~~~~~~~~~~~~~~~~ %

SamplingRate = 200;
Freqs = {Delta,Theta,Alpha,Beta,Gamma};
for f = 1:length(Freqs)
    for pat = 1:length(ProcessedHDeegFilePaths)
        fprintf('Doing patient %s...\n',num2str(pat))
        
        % ===== hdEEG =====
        fprintf('Loading hdEEG...\n')
        load(ProcessedHDeegFilePaths{pat})
        EEGnoNaN = EEGfull(:,any(~isnan(EEGfull)));
        % Narrow-bandpass filter:
        FreqL = Freqs{f}(1);
        FreqH = Freqs{f}(end);
        fprintf('Filtering...\n')
        FiltEEG = ft_preproc_bandpassfilter(bsxfun(@minus,EEGnoNaN,mean(EEGnoNaN,2)),SamplingRate,[FreqL FreqH],4,'but','twopass');
        fprintf('Calculating Hilbert transform...\n')
        Hilb = nan(size(FiltEEG));
        for c = 1:size(FiltEEG,1)
            prc_for_loop(c,size(FiltEEG,1),1);
            Hilb(c,:) = abs(hilbert(FiltEEG(c,:)));
        end
        OutputFilename = spm_file(ProcessedHDeegFilePaths{pat},'suffix',['_hilbert_',num2str(FreqL),'-',num2str(FreqH),'Hz']);
        fprintf('Saving hdEEG Hilbert...\n')
        save(OutputFilename,'Hilb','FreqL','FreqH','-v7.3')
        
        % ===== icEEG =====
        fprintf('Loading icEEG...\n')
        load(ProcessedICeegFilePaths{pat})
        % ====== Make bipolar montage (without bad channels) ======
        SEEGchannels = icChannelNames(cellfun(@isempty,regexp(icChannelNames,'Mkr')));
        SEEGchannels = SEEGchannels(cellfun(@isempty,regexp(SEEGchannels,'ecg')));
        icEEGnoNaN = icEEGd(:,any(~isnan(icEEGd)));
        [~,~,bip_labels,bip_labels2,bipSEEG] = bipolar_montage(SEEGchannels,2,icEEGnoNaN);
        fprintf('Filtering...\n')
        FiltEEG = ft_preproc_bandpassfilter(bsxfun(@minus,bipSEEG,mean(bipSEEG,2)),SamplingRate,[FreqL FreqH],4,'but','twopass');
        fprintf('Calculating Hilbert transform...\n')
        Hilb = nan(size(FiltEEG));
        for c = 1:size(FiltEEG,1)
            prc_for_loop(c,size(FiltEEG,1),1);
            Hilb(c,:) = abs(hilbert(FiltEEG(c,:)));
        end
        OutputFilename = spm_file(ProcessedICeegFilePaths{pat},'suffix',['_hilbert_',num2str(FreqL),'-',num2str(FreqH),'Hz']);
        fprintf('Saving icEEG Hilbert...\n')
        save(OutputFilename,'Hilb','FreqL','FreqH','-v7.3')
        
    end
end

% Extract epochs from Hilbert transformed signals:

for pat = 1:length(ProcessedHDeegFilePaths)
    fprintf('Doing patient %s...\n',num2str(pat))
    fprintf('Loading hdEEG...\n')
    load(ProcessedHDeegFilePaths{pat})
    load(ArtefactsInfoPaths{pat})
    clear Epochs
    
    fprintf('Loading icEEG...\n')
    load(ProcessedICeegFilePaths{pat})
    SEEGchannels = icChannelNames(cellfun(@isempty,regexp(icChannelNames,'Mkr')));
    SEEGchannels = SEEGchannels(cellfun(@isempty,regexp(SEEGchannels,'ecg')));
    [~,~,bip_labels,bip_labels2] = bipolar_montage(SEEGchannels,2);
    
    for spk = 1:length(EpochingFilePaths{pat})
        fprintf('Doing condition %s...\n',num2str(spk))
        load(EpochingFilePaths{pat}{spk})
        
        for f = 1:length(Freqs)
            FreqL = Freqs{f}(1);
            FreqH = Freqs{f}(end);
            
            load(spm_file(ProcessedHDeegFilePaths{pat},'suffix',['_hilbert_',num2str(FreqL),'-',num2str(FreqH),'Hz']))
            
            SamplingRate = NewFreq;
            FileBasename = spm_file(EpochingFilePaths{pat}{spk},'basename');
            A = nan(size(EEGfull));
            A(:,any(~isnan(EEGfull))) = Hilb;
            if strcmpi(FileBasename(end-12:end),'_spike_epochs')
                hdEpochs = get_EEG_epochs(A,SpkStart,Pre,Post);
                EpochStart = SpkStart-Pre;
                EpochEnd = SpkStart+Post;
                EpochTiming = SpikesT;
            elseif strcmpi(FileBasename(end-12:end),'random_epochs')
                hdEpochs = get_EEG_epochs(A,RandomStart,Pre,Post);
                EpochStart = RandomStart-Pre;
                EpochEnd = RandomStart+Post;
                EpochTiming = RandomsT;
            end
            
            load(spm_file(ProcessedICeegFilePaths{pat},'suffix',['_hilbert_',num2str(FreqL),'-',num2str(FreqH),'Hz']))
            
            A = nan(size(bip_labels,2),size(icEEGd,2));
            A(:,any(~isnan(icEEGd))) = Hilb;
            if strcmpi(FileBasename(end-12:end),'_spike_epochs')
                icEpochs = get_EEG_epochs(A,SpkStart,Pre,Post);
                EpochStart = SpkStart-Pre;
                EpochEnd = SpkStart+Post;
                EpochTiming = SpikesT;
            elseif strcmpi(FileBasename(end-12:end),'random_epochs')
                icEpochs = get_EEG_epochs(A,RandomStart,Pre,Post);
                EpochStart = RandomStart-Pre;
                EpochEnd = RandomStart+Post;
                EpochTiming = RandomsT;
            end
            
            save(regexprep(spm_file(EpochingFilePaths{pat}{spk},'suffix',['_hilbert_',num2str(FreqL),'-',num2str(FreqH),'Hz']),'hdEEG_','hdEEG_icEEG_'),'icEpochs','hdEpochs','FreqL','FreqH')
            
        end
        
    end
    
end
