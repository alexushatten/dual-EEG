function [ SourcesROIs, OptimalReg, SourcesSVD] = compute_inverse_in_ROIs( EEG, ISfile, ROIsFile )
%compute_inverse_in_ROIs: compute results of inverse solutions given EEG traces
% and inverse solution matrix (.is file from Cartool) and projects them to
% ROIs based on SVD method. The optimal regularization is based on L-corner
% of the Tikhonov regularization values against the average norm of all
% solution points.
%
% [SourcesROIs,OptimalReg,Sources] = compute_inverse_in_ROIs(EEG,ISfile,ROIsFile)
%
%  Inputs
% --------
% EEG: [channels x time] EEG traces
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
% See also COMPUTE_INVERSE, ESI2ROIS
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, November 2019
%-------------------------------------------------------------------------

[~,~,SourcesSVD,OptimalReg] = compute_inverse(EEG,ISfile,'optimal');
SourcesROIs = ESI2ROIs(SourcesSVD,ROIsFile);

end