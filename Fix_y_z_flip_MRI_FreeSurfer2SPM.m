
% fix strange confusion between Y and Z axis from FreeSurfer to SPM...
% Never happened before, maybe due to SPM version (SPM12b vs SPM12)...?

V = spm_vol('E:\dualEEG\sub-11\mri\T1w_class-GM_ribbon_hippo_amyg.img');
D = spm_read_vols(V);
Dout = nan(size(D));
for s = 1:size(D,1)
    Dout(s,:,:) = fliplr(squeeze(D(s,:,:))');
end
Vout = V;
Vout.fname = spm_file(Vout.fname,'prefix','rot_');
spm_write_vol(Vout,Dout);

flip_mri_left_right('E:\dualEEG\sub-11\mri\rot_T1w_class-GM_ribbon_hippo_amyg.img');

