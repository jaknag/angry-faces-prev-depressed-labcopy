# angry-faces-prev-depressed

Repository containing the input files and scripts required for data analysis for the project investigating differential processing of angry face expressions in depressed vs non-depressed people. The output files are for reference, as they are generated using the files in the input and script folders. 

input: 
- beh: folder containing .txt files which are a record of the behavioural responses made by each participant
- demog: folder containing demographics for each of the participants
-- FaceScores.csv contains the Benton and Catell face scores for each participant, arranged by CCID
-- demog_genderDiscrimination.mat contains age, gender and various other demographic parameters for each participant
-- gender_discrim_timeSinceDepression.csv contains the time (in years) since depression diagnosis for depressed participants
-- previously_depressed.csv contains a list of CCIDs of participants which have a prior history of depression
- motion: folder containing .txt files which are a record of the motion parameters
- matrices: folder containing accuracy_model_matrix1.csv which is a matrix of model parameters for automatic HDDM model estimation
- sdartel_GM_mask.nii which is the dartel grey matter mask from the task

The scripts make reference to inputs located in the following folders, which are not part of the input folder and can be obtained directly from the Cam-CAN repository:
- scans: folder containing .nii files which are the fMRI files for each participant

script: 
- script_01dataextract_outliers.m - this script extracts the behavioural data from 'beh' .txt files and excludes outliers
- script_02demographics.m - this script analyses the demographics of the participants and prepares files for further analysis
- script_03workingdemog.ipynb - this script further analyses demographics and prepares for further analysis in python
- script_04hddm.ipynb - this script fits and runs the HDDM models to the data
- script_05hddmanalysis.ipynb - this script analyses the HDDM models fitted by the previous script
- script_06importingfmri.m - this script imports and prepares the 'scans' .nii files
- script_07spm1stLevel.m - this script runs the first-level fMRI analyses using spm
- script_08spm2ndLevel.m - this script runs the second-level fMRI analyses using spm
- script_09spm2ndLevel_supp.m - this script runs the second-level fMRI analyses for supplementary material
- script_10linearmodel.R - this script runs the linear model
- script_11plotddm.R - this script plots the drift diffusion process for generation of one of the figures for the paper 

