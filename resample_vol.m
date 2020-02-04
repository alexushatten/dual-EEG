function resample_vol( ImPath, voxsiz )
% resample_vol: resample the 3D volume with specified voxel size
%
% resample_vol( ImPath, voxsiz )
%
%  Inputs
% --------
% ImPath: string, path of volume to resample
% voxsiz: [1 x 3] vector specifying new voxel size for image to resample
%
%  Outputs
% ---------
% Resampled 3D volume with voxel size specified by voxsiz.
% 'r_"v1"x"v2"x"v3"_' will be prepended to filename, with "v1", "v2" and
% "v3" being the voxel size as requested in voxsiz, respectively.
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, May 2019
%-------------------------------------------------------------------------

if length(voxsiz)~=3
    error('Specify voxsiz as a [1 x 3] vector')
end
voxsiz = voxsiz(:)'; % prevents entering [3 x 1] instead of [1 x 3]
% voxsiz = [2 2 2]; % new voxel size {mm}
% V = spm_select([1 Inf],'image');
V = spm_vol(ImPath);
for i=1:numel(V)
   bb        = spm_get_bbox(V(i));
   VV(1:2)   = V(i);
   VV(1).mat = spm_matrix([bb(1,:) 0 0 0 voxsiz])*spm_matrix([-1 -1 -1]);
   VV(1).dim = ceil(VV(1).mat \ [bb(2,:) 1]' - 0.1)';
   VV(1).dim = VV(1).dim(1:3);
   spm_reslice(VV,struct('mean',false,'which',1,'interp',0,'suffix',['_',regexprep(num2str(voxsiz),' *','x')])); % 1 for linear
end

end

