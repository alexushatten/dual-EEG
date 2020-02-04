function filtered_EEG = dual_notch_filt( EEG, Fs, Freq )
% dual_notch_filt: IIR notch filter (two-pass, 2nd order) for power line interference removal
%
% filtered_EEG = dual_notch_filt( EEG, Fs, Freq )
%
%  Inputs
% --------
% EEG  : channels x time frames EEG traces
% Fs   : sampling rate
% Freq : power line noise
% LFCO : low-frequency cut-off
%
%  Outputs
% ---------
% filtered_EEG
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, July 2018
%-------------------------------------------------------------------------

EEGmean = mean(EEG,2);
fprintf('Applying notch filter...\n');
% filtered_EEG = brainstorm_notch(bsxfun(@minus,EEG,EEGmean),Fs,Freq:Freq:LFCO);
filtered_EEG = brainstorm_notch(bsxfun(@minus,EEG,EEGmean),Fs,Freq:Freq:Fs/2); % because Butterworth does not eliminate everything 
% with the above, we eliminate also harmonics of line noise up to Nyquist
% frequency
filtered_EEG = bsxfun(@plus,filtered_EEG,EEGmean);

end

