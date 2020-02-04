function baseline_corrected_EEG = dual_filt_baseline_correct( EEG, Trials, Fs, FreqBW, FreqNotch )
% dual_filt_baseline_correct: bandpass & notch filtering, baseline correction
%
% baseline_corrected_EEG = dual_filt_baseline_correct( EEG, Trials, Fs, FreqBW, FreqNotch )
%
%  Inputs
% --------
% EEG       : channels x time frames EEG traces
% Trials    : trial timings (in samples)
% Fs        : sampling rate
% FreqBW    : frequency cut-off for Butterworth bandpass filtering
% FreqNotch : power line noise (with harmonics)
%
%  Outputs
% ---------
% baseline_corrected_EEG
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2018
%-------------------------------------------------------------------------

% % [DEFAULT VALUES EXAMPLE]
% Fs = 1000;
% NotBaseline = 25;
% FreqBW = [1 100];
% FreqNotch = [50 100 150];

filtered_EEG = dual_bw_bp_filt( EEG, Fs, FreqBW );
filtered_EEG = bsxfun(@minus,filtered_EEG,nanmedian(EEG,2)); % just to make sure no dc component left
filtered_EEG = dual_notch_filt( filtered_EEG, Fs, FreqNotch );
baseline_corrected_EEG = dual_baseline_correction( filtered_EEG, Trials, PreStim, PostStim, NotBaseline );

end

