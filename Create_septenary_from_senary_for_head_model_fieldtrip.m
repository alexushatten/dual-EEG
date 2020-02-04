
% make septenary images from senary images and region growing 3d outputs
% capturing outside-of-skull air for head model in FieldTrip:

SenaryImages = {'E:\FS_subjects_DONE\sub-01\mars\LRflp_senary.nii'};
MaskImages = {'E:\FS_subjects_DONE\sub-01\mars\LRflp_regiongrowing3d_output.nii'};

size(SenaryImages,1)
size(MaskImages,1)

for f = 1:size(SenaryImages,1)
    Vsen = spm_vol(SenaryImages{f}); Dsen = spm_read_vols(Vsen);
    Vmask = spm_vol(MaskImages{f}); Dmask = spm_read_vols(Vmask);
    Dout = Dsen.*(~Dmask);
    Vout = Vsen; Vout.fname = spm_file(Vout.fname,'basename','septenary');
    spm_write_vol(Vout,Dout);
end

