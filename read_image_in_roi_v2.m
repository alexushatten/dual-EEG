function C = read_image_in_roi_v2(Image,Mask,Idx)
% read_image_in_roi_v2: Read image intensities and coordinates in a given
% ROI (the two images should be registered but do not need same
% resolution!)
%
% S = read_image_in_roi_v2(Image,Mask,Idx)
%
%  Inputs
% --------
% Image: image file to read in
% Mask:  mask image (binary or with ROI indices)
% Idx [optional]: index to lookup if mask with multiple ROIs
%
%  Outputs
% ---------
% Structure containing image intensity at voxels in ROI, voxel and MNI
% coordinates
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, March 2019
%-------------------------------------------------------------------------

if nargin<3
    Idx = 1;
end

Mask=spm_vol(Mask);
[data_mask, XYZ]=spm_read_vols(Mask);
InMask=find(data_mask==Idx); % find voxels with intensity == 1 in the mask
XYZ=XYZ(:,InMask); % find coordinates of voxels (== 1) in the mask
XYZ(4,:)=1;
Image = spm_vol(Image);
C.XYZ = Image(1).mat\XYZ;

% % remove duplicate coordinates:
C.XYZ = str2num(char(unique(cellstr(num2str(round(C.XYZ(1:3,:)'))))))'; %#ok<ST2NM>

C.intensity=spm_get_data(Image,C.XYZ);
C.data_mask = data_mask;
C.InMask = InMask;
C.spec_mask = Mask;
C.XYZ_MNI = XYZ;
C.XYZ_MNI_final = Image(1).mat*[C.XYZ;ones(1,size(C.XYZ,2))];
C.XYZ_MNI_final = C.XYZ_MNI_final(1:3,:);

end