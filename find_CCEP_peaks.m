
load('path_to_STIMmkr_CCEP_stim_pulses_patient9.mat')

% ==== Cartool says that tracks 1, 2, 64, 93, 155 and 200 are bad when analysing the 1st stimulation run ====
% Interestingly, this is exactly what we were looking for, as 1 & 2 are
% stimulated, 64 was broken (the reference was broken and the latter
% was chosen to replace it), and the remaining ones (93, 155 and 200)
% are constantly bad throughout all recordings. The tricky thing then
% is that tracks 1 and 2 shall be removed only for the first run,
% ...etc.
BadChannels = [64, 93, 155, 200];
ECGchannels = [202,203];

PreStim = 200;
% BUT BASELINE SHOULD BE PreStim-25
NotBaseline = 25;
PostStim = 799;
HeightThreshLoose = 4; % 6 <=> 99.99966% ; 4 <=> 99.997% ; 3 <=> 99.87% (in terms of sigma)
HeightThreshStrict = 8;
HeightThresh = 6;
WidthTresh = 10; % 10 samples (equivalent to ms here)

pThreshLoose = 1/100;
pThreshStrict = 1/1000;
Ntrials = 20;

Fs = 1000;
Fco = [0.1 100]; % not appropriate for HFOs but here we focus on detection of CCEP, will adapt later for other stuff
Fco = [0.1 90]; % to match high-pass filter of Trebaul's correction method

ArtifactPeriodDuration = 10;
ArtifactPeriodDuration = 6; % for SEEG, according to Prime et al, 2017 (https://doi.org/10.1111/epi.13939)
ArtifactPeriodDuration = 0;

% Fn = Fs/2;
% Wn = Fco/Fn;
% [b1,a1] = butter(2,Wn); % 2nd order, as in Cartool
Notch = 50; % line noise (50Hz)
% Wo = Notch/Fn;
% Bw = Wo/35; % Q-factor set to 35, used to specify filter bandwidth
% [b2,a2] = iirnotch(Wo,Bw);

%% MONOPOLAR - AVG REF - corrected using optimized Trebaul's method

% New root directory:
NewRoot = 'E:\patient9_5mA_dual_custom_correction_hp90\';

UseCorrected = 1; % 0 for uncorrected signal, 1 for corrected (filtered) signal, 0.5 for corrected but unfiltered signal

[peak_lat, peak_amp, onset_lat, onset_amp, ...
    risphas_lat, risphas_amp, maxslopris_lat, maxslopris_amp,...
    offset_lat, offset_amp, desphas_lat, desphas_amp, ...
    maxslopdes_lat, maxslopdes_amp,...
    min_peak_prom_ris,min_peak_prom_des] = deal(cell(201,1,181)); % because there will be 201 channels (203 with ECG), there are 181 stimulation runs, and even if there will be more than 1 peak detected, this will help gather the results at the end of the loop

for N = 1:size(fp,1)
    MatFilePath = strrep(strrep(fp(N,:),'Z:\dualEEG\2_Resampled\patient9\stimulation\',NewRoot),'STIM_MKR_chan.mat','CCEP_5mA.mat');
    fprintf('\nStimulation run #%d, loading data...\n',N)
    load(strrep(fp(N,:),'STIM_MKR_chan.mat',['icEEG',filesep,'icEEG.mat']),'labels')
    load(MatFilePath)
    
    if UseCorrected==1
        cEEG = corrected_signal;
    elseif UseCorrected==0
        cEEG = EEG;
    elseif UseCorrected==0.5
        cEEG = corrected_unfiltered_signal;
    else
        error('Please specify whether to use corrected or uncorrected EEG')
    end
    
    %======== cleaning ========
    % clean EEG by excluding stimulated (..., ...), bad (64, 93, 155 & 200) and ECG (202 & 203) tracks:
    %     [Bipoles,ElecMatch,labelsB] = bipolar(labels);
    %     StimPairs = [Bipoles,Bipoles+1];
    % => no, the order does not match...
    % process file path instead, because it contains useful information as
    % it was systematically organized...
    [~,StimRun] = fileparts(fileparts(MatFilePath));
    TempOut1 = regexp(StimRun,'_','split');
    TempOut2 = regexp(TempOut1{3},'[0-9]+','split');
    TempOut3 = regexp(TempOut1{3},'[0-9]+');
    
    TempOut4 = regexpi(TempOut1{3},'[a-z]+');
    TempOut5 = cellfun(@length,regexpi(TempOut1{3},'[a-z]+','match'))-1;
    TempOut6 = [TempOut4(1):(TempOut4(1)+TempOut5(1)),...
        TempOut4(2):(TempOut4(2)+TempOut5(2))];
    % 'cause we know that above vectors' length can be max' 2...
    TempOut7 = TempOut1{3};
    TempOut7(TempOut6)='_';
    TempOut8 = regexp(TempOut7,'_','split');
    TempOut9 = TempOut8(~cellfun(@isempty,TempOut8))';
    
    FirstNum = TempOut9{1};
    SecondNum = TempOut9{2};
    StimChan1 = find(strcmpi([TempOut2{1},FirstNum],labels));
    StimChan2 = find(strcmpi([TempOut2{2},SecondNum],labels));
    
    %======== excluding bad channels ========
    cEEG([StimChan1,StimChan2,BadChannels, ECGchannels],:) = [];
    clabels = labels;
    clabels([StimChan1,StimChan2,BadChannels, ECGchannels])=[];
    
    %======== re-referencing ========
    % re-reference the data (average referential montage for now, later we will see for
    % bipolar montage):
    rrEEG = cEEG-repmat(mean(cEEG),size(cEEG,1),1);
    
    %======== filtering ========
    % filter (6th order fieldtrip defaults, permissive 0.5 to 100 Hz):
    fEEG = ft_preproc_lowpassfilter(rrEEG,Fs,[Fco(1) Fco(2)],6,'but','twopass');
    
    % line noise filter (4th order fieldtrip defaults):
    nEEG = notchfilter(fEEG,Fs,Notch,4);
    
    NumChan = size(nEEG,1);
    %======== baseline correction ========
    % "subtracting EACH track with ITS average within a given time period" (Cartool reference guide)
    bcEEG = nan(size(nEEG,1),Ntrials,1000); % number of electrodes x 20 epochs (stimulation pulses) x 1000 datapoints (ms)
    for n = 1:NumChan
        for epoch = 1:length(Pulses(N,:))
            bcEEG(n,epoch,:) = nEEG(n,(Pulses(N,epoch)-PreStim):(Pulses(N,epoch)+PostStim))-mean(nEEG(n,(Pulses(N,epoch)-PreStim):Pulses(N,epoch)-NotBaseline));
        end
    end
    
    %======== average ========
    aEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        aEEG(n,:) = mean(squeeze(bcEEG(n,:,:)))';
    end
    %======== deviation ========
    stdEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        stdEEG(n,:) = std(squeeze(bcEEG(n,:,:)))';
    end
    %======== t-score across trials ========
    tEEG = aEEG./(stdEEG/sqrt(Ntrials));
    %     pEEG = spm_Tpdf(tEEG,Ntrials-1);
    
    %========== matrix of significant epochs ==========
    SigMatStrict = sig_mat(tEEG, HeightThreshStrict, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMat = sig_mat(tEEG, HeightThresh, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMatLoose = sig_mat(tEEG, HeightThreshLoose, WidthTresh, PreStim, PostStim, NotBaseline);
    
    %========== detect response onset, peak and offset ==========
    % peakfit.m does not work well here... the shape of the response varies
    % a lot from one site to another!
    [pklER,slER,rplER,msrlER,elER,dplER,msdlER,...
        pkER,sER,rpER,msrER,eER,dpER,msdER,...
        FinalMinPkPromB,FinalMinPkPromE] = evoked_resp_detector(tEEG,SigMat',ArtifactPeriodDuration);
    %     % when the algorithm failed to find response start while there was a
    %     % peak:
    %     length(find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1))
    %     [I,J] = find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1);
    %     figure;
    %     for c = 1:length(I)
    %         plot(tEEG(I(c),:));
    %         hold on; plot(cell2mat(pklER(I(c),J(c))),cell2mat(pkER(I(c),J(c))),'ro');
    %         force_binary_input_from_user({},'OK to continue?');
    %     end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %========== VISUALIZE PEAKS ==========
%         figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75]);
%         subplot(2,1,1);
%         AvgStd = mean(HeightThresh*std(tEEG(:,1:PreStim-NotBaseline)')); % [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         % % --> alpha looks great, especially to overlay multiple layers of significance threshold (sigma), but when clicking on lines the figure looks a bit ugly with it...
%         [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k'); % outlinebounds(Lb,Ub);
%         % AvgStd = mean(HeightThreshLoose*std(tEEG(:,1:PreStim-NotBaseline)')); [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         hold on; set(gca,'XTick',1:50:1000); set(gca,'XTickLabel',{-200:50:800});
%     %     ylabel('amplitude (uV)');
%         ylabel('t-score'); xlabel('time (samples)'); grid on;
%         NonSigTracks = (sum(cellfun(@isempty,pklER),2)==size(pklER,2));
%         plot(tEEG'); H = get(gca,'Children');
%     %     clickText(H,vertcat(clabels,{'';'';'';''}));
%         clickText(H,flipud(vertcat({'';''},clabels))); % the two last are for the zeros (black line) and the boundedline (grey patch)
%         set(H(find(flipud(vertcat([0;0],NonSigTracks)))),'linestyle',':'); %#ok<FNDSB>
%         plot(cell2mat(pklER(~cellfun(@isempty,pklER))),cell2mat(pkER(~cellfun(@isempty,pkER))),'ko','linewidth',1,'markersize',4)
%         plot(cell2mat(slER(~cellfun(@isempty,slER))),cell2mat(sER(~cellfun(@isempty,sER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(rplER(~cellfun(@isempty,rplER))),cell2mat(rpER(~cellfun(@isempty,rpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msrlER(~cellfun(@isempty,msrlER))),cell2mat(msrER(~cellfun(@isempty,msrER))),'k^','linewidth',1,'markersize',4)
%         
%         plot(cell2mat(elER(~cellfun(@isempty,elER))),cell2mat(eER(~cellfun(@isempty,eER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(dplER(~cellfun(@isempty,dplER))),cell2mat(dpER(~cellfun(@isempty,dpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msdlER(~cellfun(@isempty,msdlER))),cell2mat(msdER(~cellfun(@isempty,msdER))),'k^','linewidth',1,'markersize',4)
%      
%         %     plot(SigEEG2,'.','linestyle','none')
%         subplot(2,1,2); imagesc(1-SigMat'); colormap('gray');
%         set(gca,'YTick',1:size(tEEG,1)); set(gca,'YTickLabel',clabels);
    %=====================================
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % re-include excluded channels as NaN channels:
    zEEG = nan(size(corrected_signal,1)-2,size(tEEG,2));
    Count = 0;
    for c = 1:size(corrected_signal,1)-2
        if ~isempty(setdiff(c,unique([Isc1;Isc2;IbcX])))
            Count = Count+1;
            zEEG(c,:) = tEEG(Count,:);
            peak_lat(c,1:size(pklER,2),N) = pklER(Count,:);
            peak_amp(c,1:size(pklER,2),N) = pkER(Count,:);
            onset_lat(c,1:size(pklER,2),N) = slER(Count,:);
            onset_amp(c,1:size(pklER,2),N) = sER(Count,:);
            
            maxslopris_lat(c,1:size(pklER,2),N) = msrlER(Count,:);
            maxslopris_amp(c,1:size(pklER,2),N) = msrER(Count,:);
            risphas_lat(c,1:size(pklER,2),N) = rplER(Count,:);
            risphas_amp(c,1:size(pklER,2),N) = rpER(Count,:);
            
            maxslopdes_lat(c,1:size(pklER,2),N) = msdlER(Count,:);
            maxslopdes_amp(c,1:size(pklER,2),N) = msdER(Count,:);
            desphas_lat(c,1:size(pklER,2),N) = dplER(Count,:);
            desphas_amp(c,1:size(pklER,2),N) = dpER(Count,:);
            
            offset_lat(c,1:size(pklER,2),N) = elER(Count,:);
            offset_amp(c,1:size(pklER,2),N) = eER(Count,:);
            
            min_peak_prom_ris(c,1:size(pklER,2),N) = FinalMinPkPromB(Count,:);
            min_peak_prom_des(c,1:size(pklER,2),N) = FinalMinPkPromE(Count,:);
        end
    end
end

PeakLat = empty_cells_to_nan_array(peak_lat);
PeakAmp = empty_cells_to_nan_array(peak_amp);
OnsetLat = empty_cells_to_nan_array(onset_lat);
OnsetAmp = empty_cells_to_nan_array(onset_amp);
RisPhaseLat = empty_cells_to_nan_array(risphas_lat);
RisPhaseAmp = empty_cells_to_nan_array(risphas_amp);
MaxSlopeRisLat = empty_cells_to_nan_array(maxslopris_lat);
MaxSlopeRisAmp = empty_cells_to_nan_array(maxslopris_amp);

OffsetLat = empty_cells_to_nan_array(offset_lat);
OffsetAmp = empty_cells_to_nan_array(offset_amp);
DesPhaseLat = empty_cells_to_nan_array(desphas_lat);
DesPhaseAmp = empty_cells_to_nan_array(desphas_amp);
MaxSlopeDesLat = empty_cells_to_nan_array(maxslopdes_lat);
MaxSlopeDesAmp = empty_cells_to_nan_array(maxslopdes_amp);

MinPeakPromRis = empty_cells_to_nan_array(min_peak_prom_ris);
MinPeakPromDes = empty_cells_to_nan_array(min_peak_prom_des);

save('C:\Users\RMAQ\Documents\MATLAB\Conduction_delays_and_amplitudes_artifact_corrected.mat','*Lat','*Amp','*_lat','*_amp','ArtifactPeriodDuration','UseCorrected','PreStim','NewRoot','NotBaseline','Fco','Notch','PostStim','HeightThresh','WidthTresh');

%% MONOPOLAR - AVG REF - not corrected for pulse artifact

% New root directory:
NewRoot = 'E:\patient9_5mA_dual_custom_correction_hp90\';
UseCorrected = 0; % 0 for uncorrected signal, 1 for corrected (filtered) signal, 0.5 for corrected but unfiltered signal

[peak_lat, peak_amp, onset_lat, onset_amp, ...
    risphas_lat, risphas_amp, maxslopris_lat, maxslopris_amp...
    offset_lat, offset_amp, desphas_lat, desphas_amp, ...
    maxslopdes_lat, maxslopdes_amp,...
    min_peak_prom_ris,min_peak_prom_des] = deal(cell(201,1,181)); % because there will be 201 channels (203 with ECG), there are 181 stimulation runs, and even if there will be more than 1 peak detected, this will help gather the results at the end of the loop

for N = 1:size(fp,1)
    MatFilePath = strrep(strrep(fp(N,:),'Z:\dualEEG\2_Resampled\patient9\stimulation\',NewRoot),'STIM_MKR_chan.mat','CCEP_5mA.mat');
    fprintf('\nStimulation run #%d, loading data...\n',N)
    load(strrep(fp(N,:),'STIM_MKR_chan.mat',['icEEG',filesep,'icEEG.mat']),'labels')
    load(MatFilePath)
    
    if UseCorrected==1
        cEEG = corrected_signal;
    elseif UseCorrected==0
        cEEG = EEG;
    elseif UseCorrected==0.5
        cEEG = corrected_unfiltered_signal;
    else
        error('Please specify whether to use corrected or uncorrected EEG')
    end
    
    %======== cleaning ========
    % clean EEG by excluding stimulated (..., ...), bad (64, 93, 155 & 200) and ECG (202 & 203) tracks:
    %     [Bipoles,ElecMatch,labelsB] = bipolar(labels);
    %     StimPairs = [Bipoles,Bipoles+1];
    % => no, the order does not match...
    % process file path instead, because it contains useful information as
    % it was systematically organized...
    [~,StimRun] = fileparts(fileparts(MatFilePath));
    TempOut1 = regexp(StimRun,'_','split');
    TempOut2 = regexp(TempOut1{3},'[0-9]+','split');
    TempOut3 = regexp(TempOut1{3},'[0-9]+');
    
    TempOut4 = regexpi(TempOut1{3},'[a-z]+');
    TempOut5 = cellfun(@length,regexpi(TempOut1{3},'[a-z]+','match'))-1;
    TempOut6 = [TempOut4(1):(TempOut4(1)+TempOut5(1)),...
        TempOut4(2):(TempOut4(2)+TempOut5(2))];
    % 'cause we know that above vectors' length can be max' 2...
    TempOut7 = TempOut1{3};
    TempOut7(TempOut6)='_';
    TempOut8 = regexp(TempOut7,'_','split');
    TempOut9 = TempOut8(~cellfun(@isempty,TempOut8))';
    
    FirstNum = TempOut9{1};
    SecondNum = TempOut9{2};
    StimChan1 = find(strcmpi([TempOut2{1},FirstNum],labels));
    StimChan2 = find(strcmpi([TempOut2{2},SecondNum],labels));
    
    %======== excluding bad channels ========
    % ==== Cartool says that tracks 1, 2, 64, 93, 155 and 200 are bad when analysing the 1st stimulation run ====
    % Interestingly, this is exactly what we were looking for, as 1 & 2 are
    % stimulated, 64 was broken (the reference was broken and the latter
    % was chosen to replace it), and the remaining ones (93, 155 and 200)
    % are constantly bad throughout all recordings. The tricky thing then
    % is that tracks 1 and 2 shall be removed only for the first run,
    % ...etc.
    cEEG([StimChan1,StimChan2,BadChannels, ECGchannels],:) = [];
    clabels = labels;
    clabels([StimChan1,StimChan2,BadChannels, ECGchannels])=[];
    
    %======== re-referencing ========
    % re-reference the data (average referential montage for now, later we will see for
    % bipolar montage):
    rrEEG = cEEG-repmat(mean(cEEG),size(cEEG,1),1);
    
    %======== filtering ========
    % filter (6th order fieldtrip defaults, permissive 0.5 to 100 Hz):
    fEEG = ft_preproc_lowpassfilter(rrEEG,Fs,[Fco(1) Fco(2)],6,'but','twopass');
    
    % line noise filter (4th order fieldtrip defaults):
    nEEG = notchfilter(fEEG,Fs,Notch,4);
    
    NumChan = size(nEEG,1);
    %======== baseline correction ========
    % "subtracting EACH track with ITS average within a given time period" (Cartool reference guide)
    bcEEG = nan(size(nEEG,1),Ntrials,1000); % number of electrodes x 20 epochs (stimulation pulses) x 1000 datapoints (ms)
    for n = 1:NumChan
        for epoch = 1:length(Pulses(N,:))
            bcEEG(n,epoch,:) = nEEG(n,(Pulses(N,epoch)-PreStim):(Pulses(N,epoch)+PostStim))-mean(nEEG(n,(Pulses(N,epoch)-PreStim):Pulses(N,epoch)-NotBaseline));
        end
    end
    
    %======== average ========
    aEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        aEEG(n,:) = mean(squeeze(bcEEG(n,:,:)))';
    end
    %======== deviation ========
    stdEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        stdEEG(n,:) = std(squeeze(bcEEG(n,:,:)))';
    end
    %======== t-score across trials ========
    tEEG = aEEG./(stdEEG/sqrt(Ntrials));
    %     pEEG = spm_Tpdf(tEEG,Ntrials-1);
    
    %========== matrix of significant epochs ==========
    SigMatStrict = sig_mat(tEEG, HeightThreshStrict, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMat = sig_mat(tEEG, HeightThresh, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMatLoose = sig_mat(tEEG, HeightThreshLoose, WidthTresh, PreStim, PostStim, NotBaseline);
    
    %========== detect response onset, peak and offset ==========
    % peakfit.m does not work well here... the shape of the response varies
    % a lot from one site to another!
    [pklER,slER,rplER,msrlER,elER,dplER,msdlER,...
        pkER,sER,rpER,msrER,eER,dpER,msdER,...
        FinalMinPkPromB,FinalMinPkPromE] = evoked_resp_detector(tEEG,SigMat',ArtifactPeriodDuration);
    %     % when the algorithm failed to find response start while there was a
    %     % peak:
    %     length(find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1))
    %     [I,J] = find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1);
    %     figure;
    %     for c = 1:length(I)
    %         plot(tEEG(I(c),:));
    %         hold on; plot(cell2mat(pklER(I(c),J(c))),cell2mat(pkER(I(c),J(c))),'ro');
    %         force_binary_input_from_user({},'OK to continue?');
    %     end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %========== VISUALIZE PEAKS ==========
%         figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75]);
%         subplot(2,1,1);
%         AvgStd = mean(HeightThresh*std(tEEG(:,1:PreStim-NotBaseline)')); % [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         % % --> alpha looks great, especially to overlay multiple layers of significance threshold (sigma), but when clicking on lines the figure looks a bit ugly with it...
%         [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k'); % outlinebounds(Lb,Ub);
%         % AvgStd = mean(HeightThreshLoose*std(tEEG(:,1:PreStim-NotBaseline)')); [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         hold on; set(gca,'XTick',1:50:1000); set(gca,'XTickLabel',{-200:50:800});
%     %     ylabel('amplitude (uV)');
%         ylabel('t-score'); xlabel('time (samples)'); grid on;
%         NonSigTracks = (sum(cellfun(@isempty,pklER),2)==size(pklER,2));
%         plot(tEEG'); H = get(gca,'Children');
%     %     clickText(H,vertcat(clabels,{'';'';'';''}));
%         clickText(H,flipud(vertcat({'';''},clabels))); % the two last are for the zeros (black line) and the boundedline (grey patch)
%         set(H(find(flipud(vertcat([0;0],NonSigTracks)))),'linestyle',':'); %#ok<FNDSB>
%         plot(cell2mat(pklER(~cellfun(@isempty,pklER))),cell2mat(pkER(~cellfun(@isempty,pkER))),'ko','linewidth',1,'markersize',4)
%         plot(cell2mat(slER(~cellfun(@isempty,slER))),cell2mat(sER(~cellfun(@isempty,sER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(rplER(~cellfun(@isempty,rplER))),cell2mat(rpER(~cellfun(@isempty,rpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msrlER(~cellfun(@isempty,msrlER))),cell2mat(msrER(~cellfun(@isempty,msrER))),'k^','linewidth',1,'markersize',4)
%         
%         plot(cell2mat(elER(~cellfun(@isempty,elER))),cell2mat(eER(~cellfun(@isempty,eER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(dplER(~cellfun(@isempty,dplER))),cell2mat(dpER(~cellfun(@isempty,dpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msdlER(~cellfun(@isempty,msdlER))),cell2mat(msdER(~cellfun(@isempty,msdER))),'k^','linewidth',1,'markersize',4)
%      
%         %     plot(SigEEG2,'.','linestyle','none')
%         subplot(2,1,2); imagesc(1-SigMat'); colormap('gray');
%         set(gca,'YTick',1:size(tEEG,1)); set(gca,'YTickLabel',clabels);
    %=====================================
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % re-include excluded channels as NaN channels:
    zEEG = nan(size(corrected_signal,1)-2,size(tEEG,2));
    Count = 0;
    for c = 1:size(corrected_signal,1)-2
        if ~isempty(setdiff(c,unique([Isc1;Isc2;IbcX])))
            Count = Count+1;
            zEEG(c,:) = tEEG(Count,:);
            peak_lat(c,1:size(pklER,2),N) = pklER(Count,:);
            peak_amp(c,1:size(pklER,2),N) = pkER(Count,:);
            onset_lat(c,1:size(pklER,2),N) = slER(Count,:);
            onset_amp(c,1:size(pklER,2),N) = sER(Count,:);
            
            maxslopris_lat(c,1:size(pklER,2),N) = msrlER(Count,:);
            maxslopris_amp(c,1:size(pklER,2),N) = msrER(Count,:);
            risphas_lat(c,1:size(pklER,2),N) = rplER(Count,:);
            risphas_amp(c,1:size(pklER,2),N) = rpER(Count,:);
            
            maxslopdes_lat(c,1:size(pklER,2),N) = msdlER(Count,:);
            maxslopdes_amp(c,1:size(pklER,2),N) = msdER(Count,:);
            desphas_lat(c,1:size(pklER,2),N) = dplER(Count,:);
            desphas_amp(c,1:size(pklER,2),N) = dpER(Count,:);
            
            offset_lat(c,1:size(pklER,2),N) = elER(Count,:);
            offset_amp(c,1:size(pklER,2),N) = eER(Count,:);
            
            min_peak_prom_ris(c,1:size(pklER,2),N) = FinalMinPkPromB(Count,:);
            min_peak_prom_des(c,1:size(pklER,2),N) = FinalMinPkPromE(Count,:);
        end
    end
end

PeakLat = empty_cells_to_nan_array(peak_lat);
PeakAmp = empty_cells_to_nan_array(peak_amp);
OnsetLat = empty_cells_to_nan_array(onset_lat);
OnsetAmp = empty_cells_to_nan_array(onset_amp);
RisPhaseLat = empty_cells_to_nan_array(risphas_lat);
RisPhaseAmp = empty_cells_to_nan_array(risphas_amp);
MaxSlopeRisLat = empty_cells_to_nan_array(maxslopris_lat);
MaxSlopeRisAmp = empty_cells_to_nan_array(maxslopris_amp);

OffsetLat = empty_cells_to_nan_array(offset_lat);
OffsetAmp = empty_cells_to_nan_array(offset_amp);
DesPhaseLat = empty_cells_to_nan_array(desphas_lat);
DesPhaseAmp = empty_cells_to_nan_array(desphas_amp);
MaxSlopeDesLat = empty_cells_to_nan_array(maxslopdes_lat);
MaxSlopeDesAmp = empty_cells_to_nan_array(maxslopdes_amp);

MinPeakPromRis = empty_cells_to_nan_array(min_peak_prom_ris);
MinPeakPromDes = empty_cells_to_nan_array(min_peak_prom_des);

save('C:\Users\RMAQ\Documents\MATLAB\Conduction_delays_and_amplitudes_artifact_uncorrected.mat','*Lat','*Amp','*_lat','*_amp','ArtifactPeriodDuration','UseCorrected','PreStim','NewRoot','NotBaseline','Fco','Notch','PostStim','HeightThresh','WidthTresh');

%% BIPOLAR - corrected using optimized Trebaul's method

% New root directory:
NewRoot = 'E:\patient9_5mA_dual_custom_correction_hp90_bipolar\';
UseCorrected = 1; % 0 for uncorrected signal, 1 for corrected (filtered) signal, 0.5 for corrected but unfiltered signal

[peak_lat, peak_amp, onset_lat, onset_amp, ...
    risphas_lat, risphas_amp, maxslopris_lat, maxslopris_amp...
    offset_lat, offset_amp, desphas_lat, desphas_amp, ...
    maxslopdes_lat, maxslopdes_amp,...
    min_peak_prom_ris,min_peak_prom_des] = deal(cell(201,1,181)); % because there will be 201 channels (203 with ECG), there are 181 stimulation runs, and even if there will be more than 1 peak detected, this will help gather the results at the end of the loop

for N = 1:size(fp,1)
    MatFilePath = strrep(strrep(fp(N,:),'Z:\dualEEG\2_Resampled\patient9\stimulation\',NewRoot),'STIM_MKR_chan.mat','CCEP_5mA.mat');
    fprintf('\nStimulation run #%d, loading data...\n',N)
    load(strrep(fp(N,:),'STIM_MKR_chan.mat',['icEEG',filesep,'icEEG.mat']),'labels')
    load(MatFilePath)
    
    if UseCorrected==1
        cEEG = corrected_signal;
    elseif UseCorrected==0
        cEEG = EEG;
    elseif UseCorrected==0.5
        cEEG = corrected_unfiltered_signal;
    else
        error('Please specify whether to use corrected or uncorrected EEG')
    end
    
    %======== cleaning ========
    % clean EEG by excluding stimulated (..., ...), bad (64, 93, 155 & 200) and ECG (202 & 203) tracks:
    %     [Bipoles,ElecMatch,labelsB] = bipolar(labels);
    %     StimPairs = [Bipoles,Bipoles+1];
    % => no, the order does not match...
    % process file path instead, because it contains useful information as
    % it was systematically organized...
    [~,StimRun] = fileparts(fileparts(MatFilePath));
    TempOut1 = regexp(StimRun,'_','split');
    TempOut2 = regexp(TempOut1{3},'[0-9]+','split');
    TempOut3 = regexp(TempOut1{3},'[0-9]+');
    
    TempOut4 = regexpi(TempOut1{3},'[a-z]+');
    TempOut5 = cellfun(@length,regexpi(TempOut1{3},'[a-z]+','match'))-1;
    TempOut6 = [TempOut4(1):(TempOut4(1)+TempOut5(1)),...
        TempOut4(2):(TempOut4(2)+TempOut5(2))];
    % 'cause we know that above vectors' length can be max' 2...
    TempOut7 = TempOut1{3};
    TempOut7(TempOut6)='_';
    TempOut8 = regexp(TempOut7,'_','split');
    TempOut9 = TempOut8(~cellfun(@isempty,TempOut8))';
    
    FirstNum = TempOut9{1};
    SecondNum = TempOut9{2};
    StimChan1 = find(strcmpi([TempOut2{1},FirstNum],labels));
    StimChan2 = find(strcmpi([TempOut2{2},SecondNum],labels));
    
    %======== excluding bad channels ========
    [Bipoles,ElecMatch,labelsB,labelsB2] = bipolar(labels);
    labelsB = labelsB';
    [Isc1, Jsc1] = find(strcmpi(labelsB2,char(labels(StimChan1))));
    [Isc2, Jsc2] = find(strcmpi(labelsB2,char(labels(StimChan2))));
    IbcX = [];
    for badc = 1:length(BadChannels)
        [Ibc, Jbc] = find(strcmpi(labelsB2,char(labels(BadChannels(badc)))));
        IbcX = [IbcX; Ibc];
    end 
    cEEG(unique([Isc1;Isc2;IbcX]),:) = []; % ECGchannels were already excluded previously when making bipolar montage
    clabels = labelsB;
    clabels(unique([Isc1;Isc2;IbcX]))=[];
    
    %======== re-referencing ========
    % was done before using Trebaul's method for artifact correction! %
    % CORRECTED ON 2018-04-09, BEFORE IT WAS SUBTRACTING THE AVERAGE FROM
    % THE EEG THAT WAS ALREADY IN BIPOLAR MONTAGE!!!! #RM
    rrEEG = cEEG;
    
    %======== filtering ========
    % filter (6th order fieldtrip defaults, permissive 0.5 to 100 Hz):
    fEEG = ft_preproc_lowpassfilter(rrEEG,Fs,[Fco(1) Fco(2)],6,'but','twopass');
    
    % line noise filter (4th order fieldtrip defaults):
    nEEG = notchfilter(fEEG,Fs,Notch,4);
    
    NumChan = size(nEEG,1);
    %======== baseline correction ========
    % "subtracting EACH track with ITS average within a given time period" (Cartool reference guide)
    bcEEG = nan(size(nEEG,1),Ntrials,1000); % number of electrodes x 20 epochs (stimulation pulses) x 1000 datapoints (ms)
    for n = 1:NumChan
        for epoch = 1:length(Pulses(N,:))
            bcEEG(n,epoch,:) = nEEG(n,(Pulses(N,epoch)-PreStim):(Pulses(N,epoch)+PostStim))-mean(nEEG(n,(Pulses(N,epoch)-PreStim):Pulses(N,epoch)-NotBaseline));
        end
    end
    
    %======== average ========
    aEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        aEEG(n,:) = mean(squeeze(bcEEG(n,:,:)))';
    end
    %======== deviation ========
    stdEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        stdEEG(n,:) = std(squeeze(bcEEG(n,:,:)))';
    end
    %======== t-score across trials ========
    tEEG = aEEG./(stdEEG/sqrt(Ntrials));
    %     pEEG = spm_Tpdf(tEEG,Ntrials-1);
    
    %========== matrix of significant epochs ==========
    SigMatStrict = sig_mat(tEEG, HeightThreshStrict, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMat = sig_mat(tEEG, HeightThresh, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMatLoose = sig_mat(tEEG, HeightThreshLoose, WidthTresh, PreStim, PostStim, NotBaseline);
    
    %========== detect response onset, peak and offset ==========
    % peakfit.m does not work well here... the shape of the response varies
    % a lot from one site to another!
    [pklER,slER,rplER,msrlER,elER,dplER,msdlER,...
        pkER,sER,rpER,msrER,eER,dpER,msdER,...
        FinalMinPkPromB,FinalMinPkPromE] = evoked_resp_detector(tEEG,SigMat',ArtifactPeriodDuration);
    %     % when the algorithm failed to find response start while there was a
    %     % peak:
    %     length(find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1))
    %     [I,J] = find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1);
    %     figure;
    %     for c = 1:length(I)
    %         plot(tEEG(I(c),:));
    %         hold on; plot(cell2mat(pklER(I(c),J(c))),cell2mat(pkER(I(c),J(c))),'ro');
    %         force_binary_input_from_user({},'OK to continue?');
    %     end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %========== VISUALIZE PEAKS ==========
%         figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75]);
%         subplot(2,1,1);
%         AvgStd = mean(HeightThresh*std(tEEG(:,1:PreStim-NotBaseline)')); % [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         % % --> alpha looks great, especially to overlay multiple layers of significance threshold (sigma), but when clicking on lines the figure looks a bit ugly with it...
%         [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k'); % outlinebounds(Lb,Ub);
%         % AvgStd = mean(HeightThreshLoose*std(tEEG(:,1:PreStim-NotBaseline)')); [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         hold on; set(gca,'XTick',1:50:1000); set(gca,'XTickLabel',{-200:50:800});
%     %     ylabel('amplitude (uV)');
%         ylabel('t-score'); xlabel('time (samples)'); grid on;
%         NonSigTracks = (sum(cellfun(@isempty,pklER),2)==size(pklER,2));
%         plot(tEEG'); H = get(gca,'Children');
%     %     clickText(H,vertcat(clabels,{'';'';'';''}));
%         clickText(H,flipud(vertcat({'';''},clabels))); % the two last are for the zeros (black line) and the boundedline (grey patch)
%         set(H(find(flipud(vertcat([0;0],NonSigTracks)))),'linestyle',':'); %#ok<FNDSB>
%         plot(cell2mat(pklER(~cellfun(@isempty,pklER))),cell2mat(pkER(~cellfun(@isempty,pkER))),'ko','linewidth',1,'markersize',4)
%         plot(cell2mat(slER(~cellfun(@isempty,slER))),cell2mat(sER(~cellfun(@isempty,sER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(rplER(~cellfun(@isempty,rplER))),cell2mat(rpER(~cellfun(@isempty,rpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msrlER(~cellfun(@isempty,msrlER))),cell2mat(msrER(~cellfun(@isempty,msrER))),'k^','linewidth',1,'markersize',4)
%         
%         plot(cell2mat(elER(~cellfun(@isempty,elER))),cell2mat(eER(~cellfun(@isempty,eER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(dplER(~cellfun(@isempty,dplER))),cell2mat(dpER(~cellfun(@isempty,dpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msdlER(~cellfun(@isempty,msdlER))),cell2mat(msdER(~cellfun(@isempty,msdER))),'k^','linewidth',1,'markersize',4)
%      
%         %     plot(SigEEG2,'.','linestyle','none')
%         subplot(2,1,2); imagesc(1-SigMat'); colormap('gray');
%         set(gca,'YTick',1:size(tEEG,1)); set(gca,'YTickLabel',clabels);
    %=====================================
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % re-include excluded channels as NaN channels:
    zEEG = nan(size(corrected_signal,1)-2,size(tEEG,2));
    Count = 0;
    for c = 1:size(corrected_signal,1)-2
        if ~isempty(setdiff(c,unique([Isc1;Isc2;IbcX])))
            Count = Count+1;
            zEEG(c,:) = tEEG(Count,:);
            peak_lat(c,1:size(pklER,2),N) = pklER(Count,:);
            peak_amp(c,1:size(pklER,2),N) = pkER(Count,:);
            onset_lat(c,1:size(pklER,2),N) = slER(Count,:);
            onset_amp(c,1:size(pklER,2),N) = sER(Count,:);
            
            maxslopris_lat(c,1:size(pklER,2),N) = msrlER(Count,:);
            maxslopris_amp(c,1:size(pklER,2),N) = msrER(Count,:);
            risphas_lat(c,1:size(pklER,2),N) = rplER(Count,:);
            risphas_amp(c,1:size(pklER,2),N) = rpER(Count,:);
            
            maxslopdes_lat(c,1:size(pklER,2),N) = msdlER(Count,:);
            maxslopdes_amp(c,1:size(pklER,2),N) = msdER(Count,:);
            desphas_lat(c,1:size(pklER,2),N) = dplER(Count,:);
            desphas_amp(c,1:size(pklER,2),N) = dpER(Count,:);
            
            offset_lat(c,1:size(pklER,2),N) = elER(Count,:);
            offset_amp(c,1:size(pklER,2),N) = eER(Count,:);
            
            min_peak_prom_ris(c,1:size(pklER,2),N) = FinalMinPkPromB(Count,:);
            min_peak_prom_des(c,1:size(pklER,2),N) = FinalMinPkPromE(Count,:);
        end
    end
end

PeakLat = empty_cells_to_nan_array(peak_lat);
PeakAmp = empty_cells_to_nan_array(peak_amp);
OnsetLat = empty_cells_to_nan_array(onset_lat);
OnsetAmp = empty_cells_to_nan_array(onset_amp);
RisPhaseLat = empty_cells_to_nan_array(risphas_lat);
RisPhaseAmp = empty_cells_to_nan_array(risphas_amp);
MaxSlopeRisLat = empty_cells_to_nan_array(maxslopris_lat);
MaxSlopeRisAmp = empty_cells_to_nan_array(maxslopris_amp);

OffsetLat = empty_cells_to_nan_array(offset_lat);
OffsetAmp = empty_cells_to_nan_array(offset_amp);
DesPhaseLat = empty_cells_to_nan_array(desphas_lat);
DesPhaseAmp = empty_cells_to_nan_array(desphas_amp);
MaxSlopeDesLat = empty_cells_to_nan_array(maxslopdes_lat);
MaxSlopeDesAmp = empty_cells_to_nan_array(maxslopdes_amp);

MinPeakPromRis = empty_cells_to_nan_array(min_peak_prom_ris);
MinPeakPromDes = empty_cells_to_nan_array(min_peak_prom_des);

save('C:\Users\RMAQ\Documents\MATLAB\Conduction_delays_and_amplitudes_artifact_corrected_bipolar.mat','*Lat','*Amp','*_lat','*_amp','ArtifactPeriodDuration','UseCorrected','PreStim','NewRoot','NotBaseline','Fco','Notch','PostStim','HeightThresh','WidthTresh');

%% BIPOLAR - not corrected for pulse artifact

% New root directory:
NewRoot = 'E:\patient9_5mA_dual_custom_correction_hp90_bipolar\';
UseCorrected = 0; % 0 for uncorrected signal, 1 for corrected (filtered) signal, 0.5 for corrected but unfiltered signal

[peak_lat, peak_amp, onset_lat, onset_amp, ...
    risphas_lat, risphas_amp, maxslopris_lat, maxslopris_amp,...
    offset_lat, offset_amp, desphas_lat, desphas_amp, ...
    maxslopdes_lat, maxslopdes_amp,...
    min_peak_prom_ris,min_peak_prom_des] = deal(cell(201,1,181)); % because there will be 201 channels (203 with ECG), there are 181 stimulation runs, and even if there will be more than 1 peak detected, this will help gather the results at the end of the loop

for N = 1:size(fp,1)
    MatFilePath = strrep(strrep(fp(N,:),'Z:\dualEEG\2_Resampled\patient9\stimulation\',NewRoot),'STIM_MKR_chan.mat','CCEP_5mA.mat');
    fprintf('\nStimulation run #%d, loading data...\n',N)
    load(strrep(fp(N,:),'STIM_MKR_chan.mat',['icEEG',filesep,'icEEG.mat']),'labels')
    load(MatFilePath)
    
    if UseCorrected==1
        cEEG = corrected_signal;
    elseif UseCorrected==0
        cEEG = EEG;
    elseif UseCorrected==0.5
        cEEG = corrected_unfiltered_signal;
    else
        error('Please specify whether to use corrected or uncorrected EEG')
    end
    
    %======== cleaning ========
    % clean EEG by excluding stimulated (..., ...), bad (64, 93, 155 & 200) and ECG (202 & 203) tracks:
    %     [Bipoles,ElecMatch,labelsB] = bipolar(labels);
    %     StimPairs = [Bipoles,Bipoles+1];
    % => no, the order does not match...
    % process file path instead, because it contains useful information as
    % it was systematically organized...
    [~,StimRun] = fileparts(fileparts(MatFilePath));
    TempOut1 = regexp(StimRun,'_','split');
    TempOut2 = regexp(TempOut1{3},'[0-9]+','split');
    TempOut3 = regexp(TempOut1{3},'[0-9]+');
    
    TempOut4 = regexpi(TempOut1{3},'[a-z]+');
    TempOut5 = cellfun(@length,regexpi(TempOut1{3},'[a-z]+','match'))-1;
    TempOut6 = [TempOut4(1):(TempOut4(1)+TempOut5(1)),...
        TempOut4(2):(TempOut4(2)+TempOut5(2))];
    % 'cause we know that above vectors' length can be max' 2...
    TempOut7 = TempOut1{3};
    TempOut7(TempOut6)='_';
    TempOut8 = regexp(TempOut7,'_','split');
    TempOut9 = TempOut8(~cellfun(@isempty,TempOut8))';
    
    FirstNum = TempOut9{1};
    SecondNum = TempOut9{2};
    StimChan1 = find(strcmpi([TempOut2{1},FirstNum],labels));
    StimChan2 = find(strcmpi([TempOut2{2},SecondNum],labels));
    
    %======== excluding bad channels ========
    [Bipoles,ElecMatch,labelsB,labelsB2] = bipolar(labels);
    labelsB = labelsB';
    [Isc1, Jsc1] = find(strcmpi(labelsB2,char(labels(StimChan1))));
    [Isc2, Jsc2] = find(strcmpi(labelsB2,char(labels(StimChan2))));
    IbcX = [];
    for badc = 1:length(BadChannels)
        [Ibc, Jbc] = find(strcmpi(labelsB2,char(labels(BadChannels(badc)))));
        IbcX = [IbcX; Ibc];
    end 
    cEEG(unique([Isc1;Isc2;IbcX]),:) = []; % ECGchannels were already excluded previously when making bipolar montage
    clabels = labels;
    clabels(unique([Isc1;Isc2;IbcX]))=[];
    
    %======== re-referencing ========
    % was done before using Trebaul's method for artifact correction! %
    % CORRECTED ON 2018-04-09, BEFORE IT WAS SUBTRACTING THE AVERAGE FROM
    % THE EEG THAT WAS ALREADY IN BIPOLAR MONTAGE!!!! #RM
    rrEEG = cEEG;
    
    %======== filtering ========
    % filter (6th order fieldtrip defaults, permissive 0.5 to 100 Hz):
    fEEG = ft_preproc_lowpassfilter(rrEEG,Fs,[Fco(1) Fco(2)],6,'but','twopass');
    
    % line noise filter (4th order fieldtrip defaults):
    nEEG = notchfilter(fEEG,Fs,Notch,4);
    
    NumChan = size(nEEG,1);
    %======== baseline correction ========
    % "subtracting EACH track with ITS average within a given time period" (Cartool reference guide)
    bcEEG = nan(size(nEEG,1),Ntrials,1000); % number of electrodes x 20 epochs (stimulation pulses) x 1000 datapoints (ms)
    for n = 1:NumChan
        for epoch = 1:length(Pulses(N,:))
            bcEEG(n,epoch,:) = nEEG(n,(Pulses(N,epoch)-PreStim):(Pulses(N,epoch)+PostStim))-mean(nEEG(n,(Pulses(N,epoch)-PreStim):Pulses(N,epoch)-NotBaseline));
        end
    end
    
    %======== average ========
    aEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        aEEG(n,:) = mean(squeeze(bcEEG(n,:,:)))';
    end
    %======== deviation ========
    stdEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        stdEEG(n,:) = std(squeeze(bcEEG(n,:,:)))';
    end
    %======== t-score across trials ========
    tEEG = aEEG./(stdEEG/sqrt(Ntrials));
    %     pEEG = spm_Tpdf(tEEG,Ntrials-1);
    
    %========== matrix of significant epochs ==========
    SigMatStrict = sig_mat(tEEG, HeightThreshStrict, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMat = sig_mat(tEEG, HeightThresh, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMatLoose = sig_mat(tEEG, HeightThreshLoose, WidthTresh, PreStim, PostStim, NotBaseline);
    
    %========== detect response onset, peak and offset ==========
    % peakfit.m does not work well here... the shape of the response varies
    % a lot from one site to another!
    [pklER,slER,rplER,msrlER,elER,dplER,msdlER,...
        pkER,sER,rpER,msrER,eER,dpER,msdER,...
        FinalMinPkPromB,FinalMinPkPromE] = evoked_resp_detector(tEEG,SigMat',ArtifactPeriodDuration);
    %     % when the algorithm failed to find response start while there was a
    %     % peak:
    %     length(find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1))
    %     [I,J] = find((cellfun(@isempty,slER)+~cellfun(@isempty,pklER))>1);
    %     figure;
    %     for c = 1:length(I)
    %         plot(tEEG(I(c),:));
    %         hold on; plot(cell2mat(pklER(I(c),J(c))),cell2mat(pkER(I(c),J(c))),'ro');
    %         force_binary_input_from_user({},'OK to continue?');
    %     end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %========== VISUALIZE PEAKS ==========
%         figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75]);
%         subplot(2,1,1);
%         AvgStd = mean(HeightThresh*std(tEEG(:,1:PreStim-NotBaseline)')); % [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         % % --> alpha looks great, especially to overlay multiple layers of significance threshold (sigma), but when clicking on lines the figure looks a bit ugly with it...
%         [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k'); % outlinebounds(Lb,Ub);
%         % AvgStd = mean(HeightThreshLoose*std(tEEG(:,1:PreStim-NotBaseline)')); [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
%         hold on; set(gca,'XTick',1:50:1000); set(gca,'XTickLabel',{-200:50:800});
%     %     ylabel('amplitude (uV)');
%         ylabel('t-score'); xlabel('time (samples)'); grid on;
%         NonSigTracks = (sum(cellfun(@isempty,pklER),2)==size(pklER,2));
%         plot(tEEG'); H = get(gca,'Children');
%     %     clickText(H,vertcat(clabels,{'';'';'';''}));
%         clickText(H,flipud(vertcat({'';''},clabels))); % the two last are for the zeros (black line) and the boundedline (grey patch)
%         set(H(find(flipud(vertcat([0;0],NonSigTracks)))),'linestyle',':'); %#ok<FNDSB>
%         plot(cell2mat(pklER(~cellfun(@isempty,pklER))),cell2mat(pkER(~cellfun(@isempty,pkER))),'ko','linewidth',1,'markersize',4)
%         plot(cell2mat(slER(~cellfun(@isempty,slER))),cell2mat(sER(~cellfun(@isempty,sER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(rplER(~cellfun(@isempty,rplER))),cell2mat(rpER(~cellfun(@isempty,rpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msrlER(~cellfun(@isempty,msrlER))),cell2mat(msrER(~cellfun(@isempty,msrER))),'k^','linewidth',1,'markersize',4)
%         
%         plot(cell2mat(elER(~cellfun(@isempty,elER))),cell2mat(eER(~cellfun(@isempty,eER))),'kv','linewidth',1,'markersize',4)
%         plot(cell2mat(dplER(~cellfun(@isempty,dplER))),cell2mat(dpER(~cellfun(@isempty,dpER))),'kx','linewidth',1,'markersize',4)
%         plot(cell2mat(msdlER(~cellfun(@isempty,msdlER))),cell2mat(msdER(~cellfun(@isempty,msdER))),'k^','linewidth',1,'markersize',4)
%      
%         %     plot(SigEEG2,'.','linestyle','none')
%         subplot(2,1,2); imagesc(1-SigMat'); colormap('gray');
%         set(gca,'YTick',1:size(tEEG,1)); set(gca,'YTickLabel',clabels);
    %=====================================
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % re-include excluded channels as NaN channels:
    zEEG = nan(size(corrected_signal,1)-2,size(tEEG,2));
    Count = 0;
    for c = 1:size(corrected_signal,1)-2
        if ~isempty(setdiff(c,unique([Isc1;Isc2;IbcX])))
            Count = Count+1;
            zEEG(c,:) = tEEG(Count,:);
            peak_lat(c,1:size(pklER,2),N) = pklER(Count,:);
            peak_amp(c,1:size(pklER,2),N) = pkER(Count,:);
            onset_lat(c,1:size(pklER,2),N) = slER(Count,:);
            onset_amp(c,1:size(pklER,2),N) = sER(Count,:);
            
            maxslopris_lat(c,1:size(pklER,2),N) = msrlER(Count,:);
            maxslopris_amp(c,1:size(pklER,2),N) = msrER(Count,:);
            risphas_lat(c,1:size(pklER,2),N) = rplER(Count,:);
            risphas_amp(c,1:size(pklER,2),N) = rpER(Count,:);
            
            maxslopdes_lat(c,1:size(pklER,2),N) = msdlER(Count,:);
            maxslopdes_amp(c,1:size(pklER,2),N) = msdER(Count,:);
            desphas_lat(c,1:size(pklER,2),N) = dplER(Count,:);
            desphas_amp(c,1:size(pklER,2),N) = dpER(Count,:);
            
            offset_lat(c,1:size(pklER,2),N) = elER(Count,:);
            offset_amp(c,1:size(pklER,2),N) = eER(Count,:);
            
            min_peak_prom_ris(c,1:size(pklER,2),N) = FinalMinPkPromB(Count,:);
            min_peak_prom_des(c,1:size(pklER,2),N) = FinalMinPkPromE(Count,:);
        end
    end
end

PeakLat = empty_cells_to_nan_array(peak_lat);
PeakAmp = empty_cells_to_nan_array(peak_amp);
OnsetLat = empty_cells_to_nan_array(onset_lat);
OnsetAmp = empty_cells_to_nan_array(onset_amp);
RisPhaseLat = empty_cells_to_nan_array(risphas_lat);
RisPhaseAmp = empty_cells_to_nan_array(risphas_amp);
MaxSlopeRisLat = empty_cells_to_nan_array(maxslopris_lat);
MaxSlopeRisAmp = empty_cells_to_nan_array(maxslopris_amp);

OffsetLat = empty_cells_to_nan_array(offset_lat);
OffsetAmp = empty_cells_to_nan_array(offset_amp);
DesPhaseLat = empty_cells_to_nan_array(desphas_lat);
DesPhaseAmp = empty_cells_to_nan_array(desphas_amp);
MaxSlopeDesLat = empty_cells_to_nan_array(maxslopdes_lat);
MaxSlopeDesAmp = empty_cells_to_nan_array(maxslopdes_amp);

MinPeakPromRis = empty_cells_to_nan_array(min_peak_prom_ris);
MinPeakPromDes = empty_cells_to_nan_array(min_peak_prom_des);

save('C:\Users\RMAQ\Documents\MATLAB\Conduction_delays_and_amplitudes_artifact_uncorrected_bipolar.mat','*Lat','*Amp','*_lat','*_amp','ArtifactPeriodDuration','UseCorrected','PreStim','NewRoot','NotBaseline','Fco','Notch','PostStim','HeightThresh','WidthTresh');


%% visualization

for N = 1:size(fp,1)
    
    figure;
    subplot(2,8,1); imagesc(empty_cells_to_nan_array(peak_lat(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(peak_lat(:,:,N)))); colorbar;
    subplot(2,8,2); imagesc(empty_cells_to_nan_array(peak_amp(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(peak_amp(:,:,N)))); colorbar;
    subplot(2,8,3); imagesc(empty_cells_to_nan_array(onset_lat(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(onset_lat(:,:,N)))); colorbar;
    subplot(2,8,4); imagesc(empty_cells_to_nan_array(onset_amp(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(onset_amp(:,:,N)))); colorbar;
    subplot(2,8,5); imagesc(empty_cells_to_nan_array(risphas_lat(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(risphas_lat(:,:,N)))); colorbar;
    subplot(2,8,6); imagesc(empty_cells_to_nan_array(risphas_amp(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(risphas_amp(:,:,N)))); colorbar;
    subplot(2,8,7); imagesc(empty_cells_to_nan_array(maxslopris_lat(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(maxslopris_lat(:,:,N)))); colorbar;
    subplot(2,8,8); imagesc(empty_cells_to_nan_array(maxslopris_amp(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(maxslopris_amp(:,:,N)))); colorbar;
    subplot(2,8,9); imagesc(empty_cells_to_nan_array(offset_lat(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(offset_lat(:,:,N)))); colorbar;
    subplot(2,8,10); imagesc(empty_cells_to_nan_array(offset_amp(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(offset_amp(:,:,N)))); colorbar;
    subplot(2,8,11); imagesc(empty_cells_to_nan_array(desphas_lat(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(desphas_lat(:,:,N)))); colorbar;
    subplot(2,8,12); imagesc(empty_cells_to_nan_array(desphas_amp(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(desphas_amp(:,:,N)))); colorbar;
    subplot(2,8,13); imagesc(empty_cells_to_nan_array(maxslopdes_lat(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(maxslopdes_lat(:,:,N)))); colorbar;
    subplot(2,8,14); imagesc(empty_cells_to_nan_array(maxslopdes_amp(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(maxslopdes_amp(:,:,N)))); colorbar;
    subplot(2,8,15); imagesc(empty_cells_to_nan_array(min_peak_prom_ris(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(min_peak_prom_ris(:,:,N)))); colorbar;
    subplot(2,8,16); imagesc(empty_cells_to_nan_array(min_peak_prom_des(:,:,N)),'alphadata',~isnan(empty_cells_to_nan_array(min_peak_prom_des(:,:,N)))); colorbar;
    
end

PL = histc(PeakLat(~isnan(PeakLat))-200,[0:10:800]);
OL = histc(OnsetLat(~isnan(OnsetLat))-200,[0:10:800]);
MSL = histc(MaxSlopeRisLat(~isnan(MaxSlopeRisLat))-200,[0:10:800]);
RPL = histc(RisPhaseLat(~isnan(RisPhaseLat))-200,[0:10:800]);

figure; bar([0:10:800],[OL,RPL,MSL,PL])
legend('onset','50% rising phase','maximal slope','peak')
% ylim([0 5500])
ylim([0 4000])
xlim([-7 805])
grid on
xlabel('Latency [ms]')
ylabel('# evoked responses')
set(gca,'xtick',[0:10:800])
set(gca,'xticklabel',{0:10:800})
box off;

figure; bar([0:10:800],[OL,RPL,MSL,PL])
legend('onset','50% rising phase','maximal slope','peak')
ylim([0 4000])
xlim([-7 805])
grid on
xlabel('Latency [ms]')
ylabel('# evoked responses')
set(gca,'xtick',[0:50:800])
set(gca,'xticklabel',{0:50:800})
box off;

%% counting / reducing dimensionality

numel(PeakLat)
sum(vectorize(isnan(PeakLat)))
sum(vectorize(~isnan(PeakLat)))
sum(vectorize((PeakLat-200)>50))
sum(vectorize((PeakLat-200)<50))
sum(vectorize((PeakLat-200)==50))
sum(vectorize((PeakLat-200)>=50))
sum(vectorize((PeakLat-200)>500))

figure; imagesc(squeeze(sum(~isnan(PeakLat),2))); colormap('hot'); colorbar;

[c,N] = find(squeeze(sum(~isnan(PeakLat),2))==8)
[c,N] = find(squeeze(sum(~isnan(PeakLat),2))==9)

% ...

figure; imagesc(cell2mat([slER(c,:);rplER(c,:);msrlER(c,:);pklER(c,:);msdlER(c,:);dplER(c,:);elER(c,:)])-200); colorbar;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%========== VISUALIZE PEAKS ==========
figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75]);
subplot(2,1,1);
AvgStd = mean(HeightThresh*std(tEEG(:,1:PreStim-NotBaseline)')); % [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
% % --> alpha looks great, especially to overlay multiple layers of significance threshold (sigma), but when clicking on lines the figure looks a bit ugly with it...
[Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k'); % outlinebounds(Lb,Ub);
% AvgStd = mean(HeightThreshLoose*std(tEEG(:,1:PreStim-NotBaseline)')); [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
hold on; set(gca,'XTick',1:50:1000); set(gca,'XTickLabel',{-200:50:800});
%     ylabel('amplitude (uV)');
ylabel('t-score'); xlabel('time (samples)'); grid on;
NonSigTracks = (sum(cellfun(@isempty,pklER),2)==size(pklER,2));
plot(tEEG'); H = get(gca,'Children');
%     clickText(H,vertcat(clabels,{'';'';'';''}));
clickText(H,flipud(vertcat({'';''},clabels))); % the two last are for the zeros (black line) and the boundedline (grey patch)
set(H(find(flipud(vertcat([0;0],NonSigTracks)))),'linestyle',':'); %#ok<FNDSB>
plot(cell2mat(pklER(~cellfun(@isempty,pklER))),cell2mat(pkER(~cellfun(@isempty,pkER))),'ko','linewidth',1,'markersize',4)
plot(cell2mat(slER(~cellfun(@isempty,slER))),cell2mat(sER(~cellfun(@isempty,sER))),'kv','linewidth',1,'markersize',4)
plot(cell2mat(rplER(~cellfun(@isempty,rplER))),cell2mat(rpER(~cellfun(@isempty,rpER))),'kx','linewidth',1,'markersize',4)
plot(cell2mat(msrlER(~cellfun(@isempty,msrlER))),cell2mat(msrER(~cellfun(@isempty,msrER))),'k^','linewidth',1,'markersize',4)

plot(cell2mat(elER(~cellfun(@isempty,elER))),cell2mat(eER(~cellfun(@isempty,eER))),'kv','linewidth',1,'markersize',4)
plot(cell2mat(dplER(~cellfun(@isempty,dplER))),cell2mat(dpER(~cellfun(@isempty,dpER))),'kx','linewidth',1,'markersize',4)
plot(cell2mat(msdlER(~cellfun(@isempty,msdlER))),cell2mat(msdER(~cellfun(@isempty,msdER))),'k^','linewidth',1,'markersize',4)

%     plot(SigEEG2,'.','linestyle','none')
subplot(2,1,2); imagesc(1-SigMat'); colormap('gray');
set(gca,'YTick',1:size(tEEG,1)); set(gca,'YTickLabel',clabels);
%=====================================
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% c =
% 
%     27
% 
% N =
% 
%     72

cell2mat(pklER(c,:))-200

% ans =
% 
%     17    42   110   299   385   494   569   644
%     

% THIS IS OK, BUT FIRST ONSET SHOULD BE LATER, SOMETHING LIKE 11 (211)
% INSTEAD OF 5 (205)

cell2mat(slER(c,:))-200

% ans =
% 
%      5    30    58   233   310   436   555   590

%% CHECK PEAKS WITH LATENCY BELOW 10 ms

sum(vectorize((PeakLat-200)<10))

% [N2check,~] = find((reshape(PeakLat,size(PeakLat,3),size(PeakLat,1),size(PeakLat,2))-200)<10)

[c2check,N2check] = find(squeeze(any(PeakLat<210,2)));

% setdiff(1:size(fp,1),N2check)'
% N2check = sort(N2check);

figure('units','normalized','outerposition',[0.1 0.1 0.75 0.75]);
Count = 0;
for N = 1:length(N2check)
    MatFilePath = strrep(strrep(fp(N2check(N),:),'Z:\dualEEG\2_Resampled\patient9\stimulation\',NewRoot),'STIM_MKR_chan.mat','CCEP_5mA.mat');
    fprintf('\nStimulation run #%d, loading data...\n',N2check(N))
    load(strrep(fp(N2check(N),:),'STIM_MKR_chan.mat',['icEEG',filesep,'icEEG.mat']),'labels')
    load(MatFilePath)
    
    if UseCorrected==1
        cEEG = corrected_signal;
    elseif UseCorrected==0
        cEEG = EEG;
    elseif UseCorrected==0.5
        cEEG = corrected_unfiltered_signal;
    else
        error('Please specify whether to use corrected or uncorrected EEG')
    end
    
    %======== cleaning ========
    % clean EEG by excluding stimulated (..., ...), bad (64, 93, 155 & 200) and ECG (202 & 203) tracks:
    %     [Bipoles,ElecMatch,labelsB] = bipolar(labels);
    %     StimPairs = [Bipoles,Bipoles+1];
    % => no, the order does not match...
    % process file path instead, because it contains useful information as
    % it was systematically organized...
    [~,StimRun] = fileparts(fileparts(MatFilePath));
    TempOut1 = regexp(StimRun,'_','split');
    TempOut2 = regexp(TempOut1{3},'[0-9]+','split');
    TempOut3 = regexp(TempOut1{3},'[0-9]+');
    
    TempOut4 = regexpi(TempOut1{3},'[a-z]+');
    TempOut5 = cellfun(@length,regexpi(TempOut1{3},'[a-z]+','match'))-1;
    TempOut6 = [TempOut4(1):(TempOut4(1)+TempOut5(1)),...
        TempOut4(2):(TempOut4(2)+TempOut5(2))];
    % 'cause we know that above vectors' length can be max' 2...
    TempOut7 = TempOut1{3};
    TempOut7(TempOut6)='_';
    TempOut8 = regexp(TempOut7,'_','split');
    TempOut9 = TempOut8(~cellfun(@isempty,TempOut8))';
    
    FirstNum = TempOut9{1};
    SecondNum = TempOut9{2};
    StimChan1 = find(strcmpi([TempOut2{1},FirstNum],labels));
    StimChan2 = find(strcmpi([TempOut2{2},SecondNum],labels));
    
    %======== excluding bad channels ========
    % ==== Cartool says that tracks 1, 2, 64, 93, 155 and 200 are bad when analysing the 1st stimulation run ====
    % Interestingly, this is exactly what we were looking for, as 1 & 2 are
    % stimulated, 64 was broken (the reference was broken and the latter
    % was chosen to replace it), and the remaining ones (93, 155 and 200)
    % are constantly bad throughout all recordings. The tricky thing then
    % is that tracks 1 and 2 shall be removed only for the first run,
    % ...etc.
    cEEG([StimChan1,StimChan2,64, 93, 155, 200, 202, 203],:) = [];
    clabels = labels;
    clabels([StimChan1,StimChan2,64, 93, 155, 200, 202, 203])=[];
    
    %======== re-referencing ========
    % re-reference the data (average referential montage for now, later we will see for
    % bipolar montage):
    rrEEG = cEEG-repmat(mean(cEEG),size(cEEG,1),1);
    
    %======== filtering ========
    % filter (6th order fieldtrip defaults, permissive 0.5 to 100 Hz):
    fEEG = ft_preproc_lowpassfilter(rrEEG,Fs,[Fco(1) Fco(2)],6,'but','twopass');
    
    % line noise filter (4th order fieldtrip defaults):
    nEEG = notchfilter(fEEG,Fs,Notch,4);
    
    NumChan = size(nEEG,1);
    %======== baseline correction ========
    % "subtracting EACH track with ITS average within a given time period" (Cartool reference guide)
    bcEEG = nan(size(nEEG,1),Ntrials,1000); % number of electrodes x 20 epochs (stimulation pulses) x 1000 datapoints (ms)
    for n = 1:NumChan
        for epoch = 1:length(Pulses(N2check(N),:))
            bcEEG(n,epoch,:) = nEEG(n,(Pulses(N2check(N),epoch)-PreStim):(Pulses(N2check(N),epoch)+PostStim))-mean(nEEG(n,(Pulses(N2check(N),epoch)-PreStim):Pulses(N2check(N),epoch)-NotBaseline));
        end
    end
    
    %======== average ========
    aEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        aEEG(n,:) = mean(squeeze(bcEEG(n,:,:)))';
    end
    %======== deviation ========
    stdEEG = nan(size(nEEG,1),1000);
    for n = 1:NumChan
        stdEEG(n,:) = std(squeeze(bcEEG(n,:,:)))';
    end
    %======== t-score across trials ========
    tEEG = aEEG./(stdEEG/sqrt(Ntrials));
    %     pEEG = spm_Tpdf(tEEG,Ntrials-1);
    
    %========== matrix of significant epochs ==========
    SigMatStrict = sig_mat(tEEG, HeightThreshStrict, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMat = sig_mat(tEEG, HeightThresh, WidthTresh, PreStim, PostStim, NotBaseline);
    SigMatLoose = sig_mat(tEEG, HeightThreshLoose, WidthTresh, PreStim, PostStim, NotBaseline);
    
    %========== detect response onset, peak and offset ==========
    % peakfit.m does not work well here... the shape of the response varies
    % a lot from one site to another!
    [pklER,slER,rplER,msrlER,elER,dplER,msdlER,...
        pkER,sER,rpER,msrER,eER,dpER,msdER,...
        FinalMinPkPromB,FinalMinPkPromE] = evoked_resp_detector(tEEG,SigMat',ArtifactPeriodDuration);
    
%     pklER(c2check(N2check(N)),:) % does not work because number of electrodes
%     changed at this stage...
    
    min(vectorize((cell2mat(pklER(~cellfun(@isempty,pklER)))-200)))
%     [c2check,~] = find((cell2mat(pklER(~cellfun(@isempty,pklER)))-200)<10)
    [c2check,~] = find(empty_cells_to_nan_array(cellfun(@(x)x<210,pklER,'uniformoutput',0))==1)
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %========== VISUALIZE PEAKS ==========
    Count = Count+1;
    subplot(6,5,Count);
    AvgStd = mean(HeightThresh*std(tEEG(:,1:PreStim-NotBaseline)')); % [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
    % % --> alpha looks great, especially to overlay multiple layers of significance threshold (sigma), but when clicking on lines the figure looks a bit ugly with it...
    [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k'); % outlinebounds(Lb,Ub);
    % AvgStd = mean(HeightThreshLoose*std(tEEG(:,1:PreStim-NotBaseline)')); [Lb,Ub] = boundedline(1:1000,zeros(1000,1),repmat(AvgStd,1,1000),'k','alpha');
    hold on; set(gca,'XTick',1:50:1000); set(gca,'XTickLabel',{-200:50:800});
    %     ylabel('amplitude (uV)');
    ylabel('t-score'); xlabel('time (samples)'); grid on;
    NonSigTracks = (sum(cellfun(@isempty,pklER),2)==size(pklER,2));
    plot(tEEG(c2check,:)'); H = get(gca,'Children');
    %     clickText(H,vertcat(clabels,{'';'';'';''}));
    clickText(H,flipud(vertcat({'';''},clabels(c2check)))); % the two last are for the zeros (black line) and the boundedline (grey patch)
    %     set(H(find(flipud(vertcat([0;0],NonSigTracks)))),'linestyle',':'); %#ok<FNDSB>
    plot(cell2mat(pklER(c2check)),cell2mat(pkER(c2check)),'ko','linewidth',1,'markersize',4)
    plot(cell2mat(slER(c2check)),cell2mat(sER(c2check)),'kv','linewidth',1,'markersize',4)
    plot(cell2mat(rplER(c2check)),cell2mat(rpER(c2check)),'kx','linewidth',1,'markersize',4)
    plot(cell2mat(msrlER(c2check)),cell2mat(msrER(c2check)),'k^','linewidth',1,'markersize',4)
    
    plot(cell2mat(elER(c2check)),cell2mat(eER(c2check)),'kv','linewidth',1,'markersize',4)
    plot(cell2mat(dplER(c2check)),cell2mat(dpER(c2check)),'kx','linewidth',1,'markersize',4)
    plot(cell2mat(msdlER(c2check)),cell2mat(msdER(c2check)),'k^','linewidth',1,'markersize',4)
    %=====================================
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    cEEG2 = EEG;
    cEEG2([StimChan1,StimChan2,64, 93, 155, 200, 202, 203],:) = [];
    rrEEG2 = cEEG2-repmat(mean(cEEG2),size(cEEG2,1),1);
    fEEG2 = ft_preproc_lowpassfilter(rrEEG2,Fs,[Fco(1) Fco(2)],6,'but','twopass');
    nEEG2 = notchfilter(fEEG2,Fs,Notch,4);
    bcEEG2 = nan(size(nEEG2,1),Ntrials,1000); % number of electrodes x 20 epochs (stimulation pulses) x 1000 datapoints (ms)
    for n = 1:NumChan
        for epoch = 1:length(Pulses(N2check(N),:))
            bcEEG2(n,epoch,:) = nEEG2(n,(Pulses(N2check(N),epoch)-PreStim):(Pulses(N2check(N),epoch)+PostStim))-mean(nEEG2(n,(Pulses(N2check(N),epoch)-PreStim):Pulses(N2check(N),epoch)-NotBaseline));
        end
    end
    aEEG2 = nan(size(nEEG2,1),1000);
    for n = 1:NumChan
        aEEG2(n,:) = mean(squeeze(bcEEG2(n,:,:)))';
    end
    stdEEG2 = nan(size(nEEG2,1),1000);
    for n = 1:NumChan
        stdEEG2(n,:) = std(squeeze(bcEEG2(n,:,:)))';
    end
    tEEG2 = aEEG2./(stdEEG2/sqrt(Ntrials));
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %========== VISUALIZE PEAKS ==========
    hold on; plot(tEEG2(c2check,:),'r:','linewidth',2)
    %=====================================
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    xlim([190 450]); % zoom on usually problematic time frames...
    
end

