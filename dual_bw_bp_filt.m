function filtered_EEG = dual_bw_bp_filt( EEG, Fs, HFCO, LFCO )
% dual_bw_bp_filt: Butterworth bandpass filtering (two-pass, 4th order)
%
% filtered_EEG = dual_bw_bp_filt( EEG, Fs, HFCO, LFCO )
%
%  Inputs
% --------
% EEG  : channels x time frames EEG traces
% Fs   : sampling rate
% HFCO : high-frequency cut-off in Hz (e.g. 1)
% LFCO : low-frequency cut-off in Hz (e.g. 100)
%
%  Outputs
% ---------
% filtered_EEG
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2018
%-------------------------------------------------------------------------

EEGmean = mean(EEG,2);
% warning off;
% /!\ WARNING SHOULD BE KEPT BECAUSE filter_with_correction.m
% WILL WARN IF AN USUAL VERSION OF filtfilt.m IS USED (THIS HAPPENS
% SOMETIMES WHEN MULTIPLE TOOLBOXES ARE ACCESSIBLE IN THE PATH (e.g.
% EEGLAB, Brainstorm, SPM, ...)) !
fprintf('Applying Butterworth filter...\n');
filtered_EEG = ft_preproc_bandpassfilter(bsxfun(@minus,EEG,EEGmean),Fs,[HFCO LFCO],4,'but','twopass');
% warning on;
filtered_EEG = bsxfun(@plus,filtered_EEG,EEGmean);

end

