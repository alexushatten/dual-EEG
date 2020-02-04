function filtered_EEG = dual_filt( EEG, Fs, FreqBW, FreqNotch )
% dual_filt: bandpass & notch filtering
%
% filtered_EEG = dual_filt( EEG, Fs, FreqBW, FreqNotch )
%
%  Inputs
% --------
% EEG       : channels x time frames EEG traces
% Fs        : sampling rate
% FreqBW    : frequency cut-off for Butterworth bandpass filtering (e.g. [1
%             100]), specified as [HighPass LowPass]
% FreqNotch : power line noise (with harmonics)
%
%  Outputs
% ---------
% filtered_EEG
%
% ===== DEFAULT VALUES EXAMPLE =====
% Fs = 1000;
% NotBaseline = 25;
% FreqBW = [1 100];
% FreqNotch = 50; % even if only 50 is specified, it will do all harmonics
% up to Nyquist frequency!
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2018
%-------------------------------------------------------------------------

%===== 2-pass 4th-order Butterworth bandpass filtering =====
filtered_EEG = dual_bw_bp_filt( EEG, Fs, FreqBW(1), FreqBW(2) );
filtered_EEG = bsxfun(@minus,filtered_EEG,nanmedian(EEG,2)); % just to make sure no dc component left
%===== 2-pass IIR 2nd-order notch filter =====
% filtered_EEG = dual_notch_filt( filtered_EEG, Fs, FreqNotch, FreqBW(2) );
filtered_EEG = dual_notch_filt( filtered_EEG, Fs, FreqNotch ); % because in some cases the line noise frequency is higher than low-frequency cutoff
%===== remove DC offset ======
filtered_EEG = detrend(filtered_EEG')';

end

