function [EEG, Time] = downsample_CIC_EEG( EEG, OldFreq, NewFreq )
% downsample_CIC_EEG: decimate (downsample) EEG using Cascaded
% Integrator-Comb filters as in Cartool software:
% 1. low-pass filter with Nyquist frequency
% cutoff
% 2. CIC filter
% 3. CIC compensation filter (slight high-pass filter)
% 4. Decimation
%
% [EEG, Time] = downsample_CIC_EEG( EEG, Hdr, NewFreq )
%
%  Inputs
% --------
% EEG: [channel x time] EEG traces
% OldFreq: original sampling rate of EEG
% NewFreq: requested sampling rate
%
%  Outputs
% ---------
% EEG: downsampled EEG
% Time: downsampled time vector
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

% LPF to prevent aliasing
EEG = ft_preproc_lowpassfilter(EEG, OldFreq, floor(NewFreq/2), 6, 'but', 'twopass');

% Time vector
Time = (0:size(EEG,2))/OldFreq;

% CIC with compensation filter and decimation
[EEG, Time] = brainstorm_resample(EEG, Time, NewFreq, 'resample-cascade');

end

