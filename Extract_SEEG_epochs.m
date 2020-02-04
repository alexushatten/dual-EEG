
SEEGfilePaths = {'file_path_to_MAT_file_with_downsampled_filtered_SEEG_signals_for_patient1'
    'file_path_to_MAT_file_with_downsampled_filtered_SEEG_signals_for_patient2'
    '...'};

ArtefactsInfoPaths = {'file_path_to_MAT_file_with_artefacts_infos_for_patient1'
    'file_path_to_MAT_file_with_artefacts_infos_for_patient2'
    '...'};

EpochsInfo = {{'path_to_MAT_file_with_epochs_infos_for_condition1_for_patient1'
               'path_to_MAT_file_with_epochs_infos_for_condition2_for_patient1'
               '...'}
              {'path_to_MAT_file_with_epochs_infos_for_condition1_for_patient2'
               'path_to_MAT_file_with_epochs_infos_for_condition2_for_patient2'
               '...'};

Pre = 199; % duration of window before event in samples (~1 second because 200 Hz)
Post = 200; % duration of window after event in samples (1 second because 200 Hz)
for pat = 1:length(SEEGfilePaths)
    fprintf('Doing patient %s...\n',num2str(pat))
    load(SEEGfilePaths{pat})
    load(ArtefactsInfoPaths{pat})
    clear Epochs
    for spk = 1:length(EpochsInfo{pat})
        fprintf('Doing condition %s...\n',num2str(spk))
        load(EpochsInfo{pat}{spk})
        SamplingRate = NewFreq;
        FileBasename = spm_file(EpochsInfo{pat}{spk},'basename');
        if strcmpi(FileBasename(end-12:end),'_spike_epochs')
            Epochs = get_EEG_epochs(icEEGd,SpkStart,Pre,Post);
            EpochStart = SpkStart-Pre;
            EpochEnd = SpkStart+Post;
            EpochTiming = SpikesT;
        elseif strcmpi(FileBasename(end-12:end),'random_epochs')
            Epochs = get_EEG_epochs(icEEGd,RandomStart,Pre,Post);
            EpochStart = RandomStart-Pre;
            EpochEnd = RandomStart+Post;
            EpochTiming = RandomsT;
        end
        
%         ChanNames = regexprep(icChannelNames,' ','none');
        
        % Write all epochs in a single .sef file
        if strcmpi(FileBasename(end-12:end),'_spike_epochs')
            OutputFile = regexprep(regexprep(regexprep(EpochsInfo{pat}{spk},'ICA_recon_with_bad_TF_i3DS2_SIWMF_avg_ref_',''),'_spike_epochs.mat','.sef'),'hdEEG_','icEEGi_');
        elseif strcmpi(FileBasename(end-12:end),'random_epochs')
            OutputFile = regexprep(regexprep(regexprep(EpochsInfo{pat}{spk},'ICA_recon_with_bad_TF_i3DS2_SIWMF_avg_ref_',''),'_random_epochs.mat','.sef'),'hdEEG_','icEEGi_');
        else
            error('End of filename does not match spike_epochs or random_epochs')
        end
        
        write_sef(OutputFile,icEEGd(:,EpochTiming)',SamplingRate,icChannelNames);
        
        % Write also the average:
        OutputFile = regexprep(OutputFile,'_epochs_','_avg_');
        write_sef(OutputFile,mean(Epochs,3)',SamplingRate,icChannelNames);
        
        % Write individual epochs in separate .sef files
        for ep = 1:size(Epochs,3)
            if strcmpi(FileBasename(end-12:end),'_spike_epochs')
                OutputFile = regexprep(regexprep(regexprep(regexprep(EpochsInfo{pat}{spk},'ICA_recon_with_bad_TF_i3DS2_SIWMF_avg_ref_',''),'epochs_',['epoch_',sprintf('%03g',ep),'_']),'_spike_epochs.mat','.sef'),'hdEEG_','icEEGi_');
            elseif strcmpi(FileBasename(end-12:end),'random_epochs')
                OutputFile = regexprep(regexprep(regexprep(regexprep(EpochsInfo{pat}{spk},'ICA_recon_with_bad_TF_i3DS2_SIWMF_avg_ref_',''),'epochs_',['epoch_',sprintf('%03g',ep),'_']),'_random_epochs.mat','.sef'),'hdEEG_','icEEGi_');
            else
                error('End of filename does not match spike_epochs or random_epochs')
            end
            write_sef(OutputFile,Epochs(:,:,ep)',SamplingRate,icChannelNames);
        end
        
    end
end

