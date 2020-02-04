# To use easy_lausanne, requirements are FreeSurfer, FSL, python v<3, numpy, scipy, nibabel, and networkx v<2
# The README states that a B0 or BOLD image should be used for coregistration, but that is only if one wants to later extract values in the "target" volume (same resolution and orientation).
# Strangely, when using a copy of the original T1 resampled at 1x1x1mm, although the Lausanne2008 atlas will have an isotropic resolution of 1mm, it will not be aligned with the anatomical image...
# The workaround is to always use (even for subjects without fMRI / diffusion acquisition) the same B0 image file resampled at 1x1x1mm for a given subject, like this:

# Run this only once after having installed easy_lausanne and configured a conda environment with relevant packages as above:
mri_convert -c -ot nii /Volumes/DATA/sub-50_b0.nii.gz /Volumes/DATA/sub-50_b0_1mm.nii.gz
# NB: do not use -cm instead of -c, because it will not necessarily be 1x1x1mm resolution, see:
# https://surfer.nmr.mgh.harvard.edu/pub/docs/html/mri_convert.help.xml.html

# Run this for every new subject:
easy_lausanne --subject_id sub-xx --target_volume /Volumes/DATA/sub-50_b0_1mm.nii.gz --target_type diffusion --output_dir /Volumes/DATA/EasyLausanne_output_sub-xx --include500
