DIR=mfilename('fullpath'); % Get the directory where this script is stored

DIR

BIDS_DIR=strcat(fileparts(DIR),filesep,'DS002');
mkdir(BIDS_DIR) %create BIDS (DS002) directory if it does not exist

SUB='sub-02-DSI';
DCM_SUB_DIR=uigetdir(fileparts(BIDS_DIR));
BIDS_SUB_ANAT=fullfile(BIDS_DIR,SUB,'anat');
mkdir(BIDS_SUB_ANAT)
system(['C:\Users\RMAQ\Documents\programs\mricrogl\dcm2niix -f ',SUB,'_T1w -o ',BIDS_SUB_ANAT,' ',DCM_SUB_DIR])
