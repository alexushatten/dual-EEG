function [rho, pval] = corr_mri( Im1, Im2, CorrType, Value )
% Estimate correlation between two MRIs across all voxels
%
% Usage:
%-------------------------------------------------------------------------
%
% [rho, p] = corr_mri( Im1, Im2, CorrType, Value )
%
%
% Inputs:
%-------------------------------------------------------------------------
%
% Im1 & Im2: path to image 1 and 2
%
% CorrType: type of correlation coefficient requested ('pearson',
% 'spearman', ...)
%
% Outputs:
%-------------------------------------------------------------------------
%
% rho: correlation coefficient value
%
% pval : p-value
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

V1 = spm_vol(Im1); V2 = spm_vol(Im2);
D1 = spm_read_vols(V1); D2 = spm_read_vols(V2);
if nargin > 3
    D1(D1~=Value)=0;D1(D1==Value)=1;
    D2(D2~=Value)=0;D2(D2==Value)=1;
end
[rho, pval] = corr(D1(:),D2(:),'type',CorrType);

end

