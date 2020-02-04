function make_chimera_senary( Im1, Im2 )
% replace scalp of first image by scalp of
% second image to create chimeric senary image,
% if tissue classes 4 (skull) are aligned
%
% Usage:
%-------------------------------------------------------------------------
%
% make_chimera_senary( Im1, Im2 )
%
% Inputs:
%-------------------------------------------------------------------------
%
% Im1: senary image based on MRI with electrodes but without gadolinium
%
% Im2: senary image based on MRI without electrodes but with gadolinium
%
%   Segmentation output with tissue classes 1-6 should be in the same
%   directory as senary images.
%
% Outputs:
%-------------------------------------------------------------------------
%
% "chimeric_senary" and "equality" images, the latter being useful for
% post-processing, e.g. mask voxels inside head keep only changes outside
% of skull
%
%-------------------------------------------------------------------------
% See also: seg2senary.m, corr_mri.m
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

CorrType = 'pearson';

Corr = corr_mri( Im1, Im2, CorrType, 4 );

if Corr > 0.5 % good alignment
    % here we will simply assume that the senary images are located in the
    % directory of segmentation output, together with tissue classes 1 to
    % 6:
    [~,~,eMRI] = my_recursive_listfiles(fileparts(Im1),'^c');
    [~,~,gMRI] = my_recursive_listfiles(fileparts(Im2),'^c');
    
    eMRI = cellstr(eMRI); gMRI = cellstr(gMRI);
    seg4fusion = eMRI(find(~cellfun(@isempty,regexp(eMRI,'c[1|2|3|4]'))));
    seg4fusion(end+1:end+2) = gMRI(find(~cellfun(@isempty,regexp(gMRI,'c[5|6]'))));
    
    V = spm_vol(char(seg4fusion));
    D = spm_read_vols(V);
    Vout = V(1);
    Vout.fname = spm_file(Vout.fname,'basename','chimeric_senary');
    [~,M] = max(D,[],4);
    Vout.pinfo = [1 0 352]';
    spm_write_vol(Vout,M);
    
    % compute and write "equality" image as well, such that post-processing
    % can be performed if desired (e.g. mask voxels inside head which
    % changed and keep only voxel changes outside of skull)
    temp = spm_read_vols(spm_vol(Im1));
    EqualTest = M==temp;
    Vnew = Vout;
    Vnew.fname = spm_file(Vnew.fname,'basename','equality');
    
    spm_write_vol(Vnew,EqualTest);
    
else
    error('Are images really aligned? Are these really the same subjects? If yes, then change the threshold... but be careful!')
end

end

