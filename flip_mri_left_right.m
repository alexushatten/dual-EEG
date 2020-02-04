function flip_mri_left_right( FilePath )
%flip_mri_left_right: flips MRI's left and right
%
% Renaud Marquis @ FBM lab, March 2018

V = spm_vol(FilePath);
D = spm_read_vols(V);
Vnew = V;
Vnew.fname = spm_file(Vnew.fname,'prefix','LRflp_');
Dnew = nan(size(D));
for s = 1:size(D,3)
    Dnew(:,:,s) = flipud(D(:,:,s));
end
spm_write_vol(Vnew,Dnew);

end