function seg2senary(seg_folder)
%seg2senary: writes an image where the maximum probability among c1, c2,
%c3, c4, c5 and c6 is coded as 1, 2, 3, 4, 5 and 6 respectively
%
% Usage:
%-------------------------------------------------------------------------
% 
% seg2senary(seg_folder)
%
% Inputs:
%-------------------------------------------------------------------------
%
% seg_folder: string with full path to folder containing segmentation
% output (c1, c2, c3, c4, c5, c6)
%
% Outputs:
%-------------------------------------------------------------------------
%
% a "senary.nii" image written in the input folder
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2018
%-------------------------------------------------------------------------

P = spm_select('ExtFPList',seg_folder,'^c.*.nii',1);
if size(P,1)~=6,error('6 images not found'),end
V = spm_vol(P);
D = spm_read_vols(V);
Vout = V(1);
Vout.fname = spm_file(Vout.fname,'basename','senary');
[~,M] = max(D,[],4);
Vout.pinfo = [1 0 352]';
spm_write_vol(Vout,M);

end