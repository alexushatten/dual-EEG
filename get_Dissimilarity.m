function [Diss, SC, GFP] = get_Dissimilarity( EEG, varargin )
% get_Dissimilarity: Dissimilarity of EEG traces, i.e. square root of the
% mean of the squared differences between the potentials measured at all
% corresponding electrodes, after Lehmann & Skrandies (1980).
%
% [Diss, GFP] = get_Dissimilarity( EEG, varargin )
%
%  Inputs
% --------
% EEG: [channels x time] EEG traces
% GFP (optional): [1 x time] GFP of EEG traces
%
%  Outputs
% ---------
% Diss: Dissimilarity of EEG between time point 1 and time point 2, etc.
%       First value of Diss is set to 0. It ranges from 0 (maps are similar)
%       to 2 (maps are inverted)
% SC: "Spatial" correlation, defined after Brandeis et al. (1992)
% GFP: GFP of EEG traces
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

if nargin > 1
    GFP = varargin{1};
else
    [GFP,AvgRefEEG] = get_GFP(EEG);
end

% % re-reference:
% EEG = bsxfun(@minus,EEG,mean(EEG,1));

% shift by 1 sample:
EEGoverGFP = bsxfun(@rdivide,AvgRefEEG,GFP);
% U = bsxfun(@rdivide,EEG(:,1:end-1),GFP(1:end-1));
% V = bsxfun(@rdivide,EEG(:,2:end),GFP(2:end));

% calculate dissimilarity:
Diss = zeros(size(GFP));
Diss(2:end) = sqrt(sum((EEGoverGFP(:,2:end)-EEGoverGFP(:,1:end-1)).^2)/size(EEG,1));
% Diss(2:end) = sqrt(sum((V-U).^2)/size(EEG,1));

if nargout > 1
    SC = (2-Diss.^2)/2;
end

end

