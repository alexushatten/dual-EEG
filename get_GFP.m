function [GFP, AvgRefEEG] = get_GFP( EEG )
% get_GFP: Global Field Power of EEG traces
%
% [GFP, AvgRefEEG] = get_GFP( EEG )
%
%  Inputs
% --------
% EEG: [channels x time] EEG traces
%
%  Outputs
% ---------
% GFP: GFP of EEG traces
% AvgRefEEG: average-referenced EEG
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

AvgRefEEG = bsxfun(@minus,EEG,mean(EEG,1));
GFP = sqrt(sum(AvgRefEEG.^2)/size(EEG,1));
% NB: var(EEG,1,1) would give the same but is slightly slower

end

