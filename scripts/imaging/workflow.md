### General workflow for neuroimaging analyses

##QC check for VIF outliers

make sure there is a /results/imaging/STmodel3_cues_vifsQC folder
run fMRI_QC.m (change the thresh value. Plot will be produced in local to visualize vifs imgs/sub) 
	saves subjectVifs.mat structure
	contains:
		.subjects names (based in func imgs folder) and vifs metrics for each contrast img
		.outvifs containing vector of 0 or 1 (where vif>thresh) in each cell/sub 
		.cont : names of contrast selected
		.vifs : vif vector for each contrast, output of scn_spm_design_check() fct
	-contrast are selected from the SPM.mat in each subject folder, and the string to extract contrast is 'AllPains_trial0'
	 to get cue evoked contrast only. Note that the begining (to extract data...) of the code is the same as in run_med.m

## Mediation
run run_med.m
	change model and variables accordingly
	-Requires 'SCEBLmri_learndata_FINAL_N36.mat' (learn data) and 'finalPreproc_simGen_data.mat' (genData), and 
