::1st argument: subject ID
::2nd argument: data type (anat, dwi, func, fmap)
::3rd argument: input directory
::4th argument: output directory
@echo off
::filepath of script with filename and extension
::SET DIR=%0

SET DCM2NIIX_DIR="C:\Program Files\mricrogl"
:: change to correct drive (adapt if dcm2nii is not on C:\)
C:

::folder containing the script:
SET DIR=%~dp0
SET DIR=%4
::echo %DIR%

SET SUBJ_DIR=%1
SET DATA_TYPE_DIR=%2
SET INPUT_DIR=%3
SET ID=%1
::SET OUTPUT_SUBFOL=%4

::SET BIDS_DIR=%DIR%NewFolder\Test
SET BIDS_DIR=%DIR%\%SUBJ_DIR%\%DATA_TYPE_DIR%
::echo %BIDS_DIR%

mkdir %BIDS_DIR%

cd %DCM2NIIX_DIR%

::Run command for conversion of DICOM to NIfTI using dcm2niix
echo ============================================= >> %DIR%\"log.txt"
echo dcm2niix -b y -z y -v y -o %BIDS_DIR% -f "%ID%_%%d_%%t_%DATA_TYPE_DIR%" %INPUT_DIR% >> %DIR%\"log.txt"
dcm2niix -b y -z y -v y -o %BIDS_DIR% -f "%ID%_%%d_%%t_%DATA_TYPE_DIR%" %INPUT_DIR% >> %DIR%\"log.txt"
type %DIR%\"log.txt"

::dcm2niix -b y -z y -v y -o %BIDS_DIR% -f "%SUBJ_DIR%_%DATA_TYPE_DIR%_%%t_%%d" %INPUT_DIR%
::dcm2niix -b y -z y -v y -o %BIDS_DIR% -f "%SUBJ_DIR%_%DATA_TYPE_DIR%" %INPUT_DIR%


::echo "dcm2niix -b y -z y -v y -o %BIDS_DIR% -f "%SUBJ_DIR%_%DATA_TYPE_DIR%_%%t_%%d" %INPUT_DIR%"
::echo "dcm2niix -b y -z y -v y -o %BIDS_DIR% -f "%SUBJ_DIR%_%DATA_TYPE_DIR%" %INPUT_DIR%"
