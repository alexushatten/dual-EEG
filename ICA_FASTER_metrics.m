function [ MeanGradientHF, SpatialKurtosis, IC_hurst_exp, Zscores ] = ICA_FASTER_metrics( Activations, Winv )
% ICA_FASTER_metrics: estimate independent components FASTER metrics, based
% on Nolan et al (2010), with some modifications to the original pipeline.
% Because we use a bandpass Butterworth filter, we don't look at the mean
% slope around low-pass filter band for each component.
%
% [ MeanGradientHF, SpatialKurtosis, IC_hurst_exp ] = ...
%                                   ICA_FASTER_metrics( Activations, Winv )
%
%  Inputs
% --------
% Winv : [channels x components] mixing matrix
% Activations : [components x time] time courses of components
%
%  Outputs
% ---------
% MeanGradientHF: 
% SpatialKurtosis: 
% IC_hurst_exp: 
% Zscores: 
%
%-------------------------------------------------------------------------
% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

% NumSamp = size(EEGd,2)/NewFreq*SamplingFreq;
% Fn = SamplingFreq/2;
% [B, A] = butter(4, [LowPassBand/Fn HighPassBand/Fn]);
% % [B, A] = butter(4, HighPassBand/Fn);
% bw = noisebw(B,A,NumSamp,SamplingFreq)
% 
% lpf_band = 

% Here we do not look at mean slope around low-pass filter band, because we
% used a Butterworth bandpass, so we will just ignore this part.

% #RM@FBMlab, NB: correlation with EOG channels is done separately, not in this
% function as in the original method.

[MeanGradientHF,SpatialKurtosis,IC_hurst_exp] = deal(nan(size(Activations,1),1));
for c = 1:size(Activations,1)
    % temporal properties:
    MeanGradientHF(c) = mean(diff(Activations(c,:)));
    % spatial properties:
    SpatialKurtosis(c) = kurt(Winv(:,c));
    % other properties:
    IC_hurst_exp(c) = hurst_exponent(Activations(c,:));
end

Zscores = normalize_by_col(MAD_normalize_by_col([MeanGradientHF,SpatialKurtosis,IC_hurst_exp]));

end

