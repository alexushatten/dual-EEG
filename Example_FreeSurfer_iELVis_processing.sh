
############## FreeSurfer preprocessing of pre-implant MRI ##############

# make sure $SUBJECTS_DIR is writable:
sudo chmod -R 777 /usr/local/freesurfer/subjects

# process MRI data using FreeSurfer through all steps of cortical reconstruction:
recon-all -s sub-36 -i /DATA/sub-36/anat/sub-36_anat.nii.gz -all -localGI
# add -openmp INT at the end to allocate ressources for multiple cores, but do not run as many subjects as number of cores -1 simultaneously, because there is a step in recon-all which likely uses hyperthreading and can freeze processing...

# this will create a "sub-36" directory in $SUBJECTS_DIR with all necessary folders and start cortical reconstruction which will be veeeeeery long...

# Sometimes, FreeSurfer will fail to process subject, this can sometimes be fixed by running the following
recon-all -s sub-36 -hemi rh -autorecon2 -autorecon3
recon-all -s sub-36 -hemi lh -autorecon2 -autorecon3
recon-all -s sub-36 -autorecon2-wm -hemi lh -randomness
recon-all -s sub-36 -localGI

# to check the accuracy of the automated segmentation:
tkmedit sub-36 T1.mgz -aseg -surfs

# to check the accuracy of the pial surface:
tksurfer sub-36 lh pial
tksurfer sub-36 rh pial

############## coregister post-implant CT to MRI ##############
ct2mri.sh sub-36 /DATA/sub-36/ct/CT.nii.gz
# if using post-implant MRI then use mri2mri.sh instead

# T1.nii.gz : full head MRI
# brainmask.nii.gz: The skull stripped MRI
# ctINt1.nii.gz: The post-implant CT coregistered to the pre-implant MRI (if you used ct2mri.sh)
# postT1INpreT1.nii.gz: The post-implant MRI coregistered to the pre-implant MRI (if you used mri2mri.sh)

# images will automatically generated to check the quality of coregistration

# run command below to interactively inspect coregistration:
fslview /usr/local/freesurfer/subjects/sub-36/elec_recon/T1.nii.gz /usr/local/freesurfer/subjects/sub-36/elec_recon/ctINt1.nii.gz

############## localize electrodes in the coregistered CT scan using BioImageSuite ##############

# Start BioImage Suite (http://bioimagesuite.yale.edu/manual/guide/bioimagesuite_manual_95522_284_16885_v1.pdf).

/usr/local/bioimagesuite30_64/start_bioimagesuite

# In the menu that pops up, click the Editors tab on the left and then select Electrode Editor in the new main menu that appears.
# Two new windows will pop up. The Electrode Editor window is for viewing neuroimaging files (your CT and MR) and the
# Electrode Control window is for creating the *.mgrid text file that indicates where electrodes are.

# In the Electrode Editor window, load the CT or MRI nii file via the File menu.
# If using a CT scan, go to the Electrode Editor window and change view from 3-slice mode
# to 3D only (upper button) & 3-slice to Volume (below and left) to see the CT scan as 3D volume instead of slices.
# At first you will just see the CT scan as a block or cylinder because it is filling in every voxel.
# If using a postimplant MRI, select a 2D view (e.g., "Axial only").

# If using a CT scan, in the Electrode Editor window,
# click on Image Processing->Threshold and a new threshold window will pop up.
# Enter something like 5000 as the low threshold and something like 15000 or 20000
# as the high threshold as default values. The best threshold value may depending on
# your CT scanner and particular data set. Now you should see the electrodes (and maybe teeth and wires).
# You may need to play with these thresholds a bit to improve the definition of electrodes (especially
# for electrodes that are almost overlapping). If electrodes look pixilated, try raising the upper limit.

# In the Electrode Control window, under the Patient Info tab, you can add or delete sets of electrodes.
# By default, it starts with an 8x8 grid. In the Electrode Control window, under the Electrode Info tab you
# can edit the properties of sets of electrodes. Probably the only two fields you need to edit are Grid Name
# and and Dimensions. When naming the set of electrodes, be sure to use the naming conventions described
# below (see CRITICAL!) so the name indicates the electrode type and hemisphere.
# Once you change any of these properties, hit the Update Grid button and indicate Yes in the pop up dialogue box.

# To be able to start assigning electrodes to the postimplant CT  or MRI, activate Button Pick and
# press Pick Current in the lower right corner of the Electrode Info tab on the Electrode Control window.
# You also need to switch to Full Edit Mode under the Edit menu of the Electrode Control window.
# Change the color of the set of electrodes to something unique and easy to see via Grid->Grid Color on
# the Electrode Control window menu.

# To assign a particular electrode to a blob in the postimplant scan, select that electrode in the Electrode Arrangement
# section of the Electrode Control window. Under Electrode Properties the Grid/Electrode line will tell you the number
# of the set of electrodes (i.e., strip or grid; indexing starts with 0) and the number of the electrode in that set
# that you've selected (indexing starts with 1). Hold down "shift" and left click on the electrode in the CT/MRI scan.
# You may need to hold down the mouse button for a bit. A colored sphere or circle should appear where you're pointing.
# If it doesn't seem to work, switch to another electrode in that grid/strip/depth and come back to it.

# If you're not happy with the localization in a 3D CT scan try changing Volume Trace to a lower value (default=10) via
# the Volume Trace menu in the Electrode Control window. This parameter specifies the size of the sphere (I believe it
# defines the diameter in mm) around the mouse pointer that is averaged over to determine the center of mass where the
# electrode is placed. Smaller values means the electrode center will be closer to where the mouse is located.
# Raising the CT threshold may also help as faint above threshold voxels may capture the electrode.
# Raising the CT scan threshold can also help as some voxels with low CT signal may not be visible to your eye,
# but will pull the electrode off center. I do not think Volume Trace has any effect when indicating electrode
# locations in 2D slices.

# When you're done identifying all the electrodes, save the locations to a text file via the File->Save menu on
# the Electrode Control window. Call it *.mgrid where * is the patient's codename (e.g., PT001.mgrid).
# Note, that BioImage Suite automatically saves your results as you working in an *autosave.mgrid file.

#======= IMPORTANT ! ======= 

# Start each electrode name with RD_, RS_, RG_, LD_, LS, or LG_ followed by the name of the electrode (e.g., RD_RHD).
# The first letter indicates which hemisphere the electrode lies in and the second letter indicates if the
# electrode is a depth, subdural strip, or subdural grid electrode. For example, RD_RHD might be a name of
# right hemisphere hippocampal depth electrode and LG_Grid might be a name of left hemisphere grid of subdural electrodes.

# Do NOT put spaces or additional underscores in the names of the electrodes

#======= USEFUL TIPS ======= 

# You can select how many sets of electrodes to view (one, all, a subset) via the Display menu of the Electrode Control window.

# It can be helpful to use the Volume button on the Electrode Editor window to selectively hide parts of the CT scan.

# You can get rid of the green box around the CT scan with the Box menu on the Electrode Editor window.

# To view electrodes on the MRI, load the MRI nii file via the File menu of the Electrode Editor window.

# The right Zoom arrow on the Electrode Editor window zooms in, the left zooms out.

# Rotate the image in the Electrode Editor window by clicking and moving your mouse on the image.

# Click the Reset button on the Electrode Editor window to return to default viewing angle and scale.

# Export the CT/MRI image to disk with the Save Snapshot button on the Electrode Editor window.

# You can also label electrode locations in 2D slice views.
# This is necessary for labeling electrodes in postimplant MRIs and useful for electrodes that show up only weakly in the CT.

# If it is difficult to label to label all the electrodes in a CT scan (due to poor CT resolution or
# small electrodes/electrode spacing), you can label 3 contacts and then select menu Grid > Auto Warp.
# BioimageSuite will then fit a spline through the points, positioning all the other electrodes along the shaft automatically.
# This has the advantage of "smoothing" the electrode as compared to manually clicking on every artifact.

# For strips and depths, the electrode furthest from the opening in the skull is #1.

# To see the numbers of the electrodes in the Electrode Editor window, select a non-zero fontsize in the Labels menu of
# the Electrode Control window.

# Try to give each electrode a unique color. Use colors that will show up well on the MRI.
# Thus, it is best to avoid gray scale colors. It is best to use bright (e.g., neon) colors for depths.

############## export mgrid information and correct for postimplant brain shift ##############

# The electrode coordinates in the mgrid file you created need to be exported to simple text files and
# converted to a few different coordinate systems for plotting. We do this with one of the following two
# sets of functions that also correct subdural electrode locations for postimplant brain shift.
# Depth electrode coordinates are NOT corrected for brain shift.

# Start Matlab and use either Yang, Wang et al (2011)'s method:

# Make *PostimpLoc.txt in patient's elec_recon folder
matlab -nosplash -nodesktop -r 'makeIniLocTxtFile('sub-36'); quit'

# Make multiple coordinate files (including coordinates corrected for brain shift)
matlab -r 'yangWangElecPjct('sub-36');'

# This method does NOT work on grids placed over medial brain areas (e.g., an interhemispheric grid).
# To use this method on such grids, simply relabel the grids as strips in the BioImage Suite mgrid file.

# ... or Dykstra et al (2011)'s method:

# Make multiple coordinate files (including coordinates corrected for brain shift)
matlab -r 'dykstraElecPjct('sub-36');'

############## create images of electrode locations ##############

matlab -r 'plotMgridOnPial('sub-36',1);'
matlab -r 'cfg=[]; cfg.printFigs=1; plotMgridOnSlices('PT001',cfg);'

############## get anatomical label index ##############
matlab -r 'parcOut=elec2Parc('sub-36','D');'
matlab -r 'parcOut=elec2Parc('sub-36','DK');'
# Yeo et al (2011) atlas also available but need to create the mapping first with:
matlab -r 'createIndivYeoMapping('PT001'); quit'

############## get Proximal Tissue Density index ##############

matlab -r 'PTD_idx = getPtdIndex('sub-36');'


