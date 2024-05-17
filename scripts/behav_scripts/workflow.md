## General workflow for pain generalization project

### Preparation scripts and variables preprocessing
Run LKscripts/analyze_main_SCEBL_fMRI_only36FinalSubjects.m 
	-Change projectDir
	output : save in /results/behavioral/SCEBLmri_data_FINAL_N36.mat

run analyze_simrate.m (~/LKscripts)
	outputs/~results/behavioral/simmat.mat  simrate.mat 
	 # modified simrate, to export simmat3 (order of cond diff)

run SCEBmri_gen_SCEBL_fMRI.m (in ~scripts/LKscripts)
	output ~/results/behavioral/SCEBLmri_data_gen36.mat (SCEBLmri_gendata variable)

Download/add 'smDistance-wingfield2022.xlsx' to 'results/behavioral/'
	-comes from Landcaster U, sensorimotor calculator, sim matrix of learn cues

Run simModels_LearnGenData_asCSV.m
	-Specify similarity models, compute trial simratings
	-saves up-to-date Trial x Trial in a .csv at /results/SCEBLmri_Finaldata_TxT_N36.csv 
	- also process trial by trial sim ratings, as well as
	-  median sim matrices for different models (e.g. Landcaster, theoretical, empirical...)
	- Save SCEBLmri_gendata.mat (with added var) to ~/results

### Main effects

run LKscriptOutput_mainEffects.m
	-wil run multilevel GLM for learning and generalization task
	- plot used for NeurAI congress (mean diff. of pain ~ CS) as well

## fMRI analysis

### QUality check and outlier removals
run fMRI_QC.m
	- Will perform 1)func imgs extraction, 2)vifs extractions 3)remove trials above thresh
	- Change thresh at begining of script accordingly
	- Will save a subject struct having subj names, vifs per trials
		 <save ('subjectsVifs.mat', "subjects")> 
	- vifStruct.outvifs contains vector of trials exceeding vif_thresh (1 if >, 0 if <)

run run_med.m 
	- performs data extraction, reorganize files in cell/subs, mean centering vars..
	- removes trials for all variables for which the vifStruct.outvifs has 1 in the vector
	- performs multilevel mediation with boot strap...
	- change saveDir, and it will save a folder/model in saveDir

run cv_pred.m


 
## Extract outputs from ~/results to perform main contrasts
run /behav_scripts/LKscriptOutput_mainEffects.m
	save nothing. Use to display and look at multilevel behavioral gen results






# Useful code

### np.unique(x, return_counts=True)
[uniqueValues, ~, idx] = unique(x);
counts = histcounts(idx, 1:numel(uniqueValues)+1);

% Display unique values and counts
disp('Unique Values:');
disp(uniqueValues');
disp('Counts:');
disp(counts');

