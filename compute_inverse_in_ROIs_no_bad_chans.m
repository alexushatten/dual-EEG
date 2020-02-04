function [ SourcesROIs, OptimalReg, SourcesSVD] = compute_inverse_in_ROIs_no_bad_chans( Data, ISfile, ROIsFile )
%compute_inverse_in_ROIs_no_bad_chans: compute results of inverse solutions given EEG traces
% and inverse solution matrix (.is file from Cartool) and projects them to
% ROIs based on SVD method. The optimal regularization is based on L-corner
% of the Tikhonov regularization values against the average norm of all
% solution points. Bad channels are ignored.
%
% [SourcesROIs,OptimalReg,Sources] = compute_inverse_in_ROIs(Data,ISfile,ROIsFile)
%
%  Inputs
% --------
% Data: structure with fields:
%                   - EEG: [channels x time] EEG traces
%                   - bad: [n x 1] vector with indices of bad channels
% ISfile: char, path to inverse solution file (.is) (output of Cartool)
% ROIsFile: char, path to region of interest file (.rois) (output of Cartool)
%
%  Outputs
% ---------
% SourcesROIs: first eigenvector of results of inverse solution in solution
%          	   ROI space
% OptimalReg: index of optimal regularization (based on L-corner)
% Sources: first eigenvector of results of inverse solution in solution
%          points space
%
% See also COMPUTE_INVERSE, ESI2ROIS, COMPUTE_INVERSE_IN_ROIS
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2019
%-------------------------------------------------------------------------

[~,~,SourcesSVD,OptimalReg] = compute_inverse_no_bad_chans(Data,ISfile,'optimal');
SourcesROIs = ESI2ROIs(SourcesSVD,ROIsFile);

end