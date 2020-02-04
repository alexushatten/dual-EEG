::1st argument: subject ID
::2rd argument: input directory
::3rd argument: output directory
@echo off
::filepath of script with filename and extension
::SET DIR=%0

SET DCM2NII_DIR="C:\Program Files\mricron"
:: change to correct drive (adapt if dcm2nii is not on C:\)
C:

::folder containing the script:
SET DIR=%~dp0
SET DIR=%3
::echo %DIR%

SET SUBJ_DIR=%1
SET DATA_TYPE_DIR=ct
SET INPUT_DIR=%2
SET ID=%1
::SET OUTPUT_SUBFOL=%4

::SET BIDS_DIR=%DIR%NewFolder\Test
SET BIDS_DIR=%DIR%\%SUBJ_DIR%\%DATA_TYPE_DIR%
::echo %BIDS_DIR%

mkdir %BIDS_DIR%

cd %DCM2NII_DIR%

:: Run command for conversion of DICOM to NIfTI using dcm2niix
echo ============================================= >> %DIR%\"log.txt"
echo dcm2nii -t y -o %BIDS_DIR% %INPUT_DIR% >> %DIR%\"log.txt"
dcm2nii -t y -o %BIDS_DIR% %INPUT_DIR% >> %DIR%\"log.txt"
type %DIR%\"log.txt"
:: copy "dcm2niigui.ini" %DIR%\"dcm2nii_settings.txt"

:: Delete reoriented image (latest output)
echo Deleting re-oriented image...
dir %BIDS_DIR% /b /o-d > C:\TEMP\temp123456.txt && for /l %%l in (1,1,1) do @for /f "tokens=1,2* delims=:" %%a in ('findstr /n /r "^" C:\TEMP\temp123456.txt ^| findstr /r "^%%l:"') do @del %BIDS_DIR%\%%b
del C:\TEMP\temp123456.txt


