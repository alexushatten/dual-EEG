
:: Template Windows batch file that compresses and converts DICOM data using dcm2niix (and fallback dcm2nii for CT scans)

::============================ COMPRESSION ============================
:: this is useful for sharing anonymized data with all DICOM flags (the only way to make sure where is left and right)

"C:\Program Files\7-Zip\7z.exe" a E:\anonymized_DICOMs.zip E:\anonymized_DICOMs
:: Do this for anatomical, diffusion and CT data
:: NB: ignore gadolinium-enhanced T1 as it poses problem, but keep it when MPRAGE sequence has been performed with EEG cap, because glueing electrodes on scalp will be better without electrode artefacts
:: NB: no good contrast between gray and white matter in INV2 for mp2rage

::============================ CONVERSION ============================

::~~~~~~~~~~~~ CONVERT T1w (MPRAGE AND/OR MP2RAGE) DICOMs ~~~~~~~~~~~~
C:\dicom2niftis.bat sub-XX anat E:\path_to_DICOMs_anatomical_MRI E:\dcm2niix_output\sub-XX

::~~~~~~~~~~~~ CONVERT DWI DICOMS ~~~~~~~~~~~~
C:\dicom2niftis.bat sub-XX dwi E:\path_to_DICOMs_diffusion_MRI E:\dcm2niix_output\sub-XX

::~~~~~~~~~~~~ CONVERT CT DICOMs ~~~~~~~~~~~~
C:\ct_dicom2niftis.bat sub-02 E:\path_to_DICOMs_CT_scan E:\dcm2niix_output\sub-XX\

