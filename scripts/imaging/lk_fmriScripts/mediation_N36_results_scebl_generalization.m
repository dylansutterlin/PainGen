%% Thresholding with FDR

clear all
close all
clc

% mask = which('brainmask.nii');
mask = which('gray_matter_mask.img');

%% [FDR corrected] mediation results similarity to CShigh during CUE

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/mediation_38/SimHC_WBgencues_painrate_36

SETUP = mediation_brain_corrected_threshold('fdr', 'mask', 'mask.img');  % to get FDR threshold across all mediation result images  

[clpos_gencue_A_001, clneg_gencue_A_001, clpos_data_gencue_A_001, clneg_data_gencue_A_001] = mediation_brain_results('a', 'thresh', .001, 'size', 5, 'mask', mask)
[clpos_gencue_A, clneg_gencue_A, clpos_data_gencue_A, clneg_data_gencue_A] = mediation_brain_results('a', 'thresh', [.001 .01], 'size', [5 20], 'mask', mask)
[clpos_gencue_B_001, clneg_gencue_B_001, clpos_data_gencue_B_001, clneg_data_gencue_B_001] = mediation_brain_results('b', 'thresh', .001, 'size', 5, 'mask', mask)
[clpos_gencue_B, clneg_gencue_B, clpos_data_gencue_B, clneg_data_gencue_B] = mediation_brain_results('b', 'thresh', [.001 .01], 'size', [5 20], 'mask', mask)
[clpos_gencue_AB_001, clneg_gencue_AB_001, clpos_data_gencue_AB_001, clneg_data_gencue_AB_001] = mediation_brain_results('ab', 'thresh', .001, 'size', 5, 'mask', mask)
[clpos_gencue_AB, clneg_gencue_AB, clpos_data_gencue_AB, clneg_data_gencue_AB] = mediation_brain_results('ab', 'thresh', [.001 .01], 'size', [5 20], 'mask', mask)

xmy_gencue = fmri_data('X-M-Y_effect.img');

%% make plots for CUE mediators 

gencue_disp = fmridisplay;
% gencue_A_disp = montage(gencue_A_disp, 'axial', 'wh_slice', [0 0 -66; 0 0 -10; 0 0 10; 0 0 30; 0 0 40; 0 0 50], 'onerow');
gencue_disp = montage(gencue_disp, 'axial', 'wh_slice', [0 0 -66; 0 0 -10; 0 0 0; 0 0 10; 0 0 20; 0 0 30; 0 0 40; 0 0 50; 0 0 60], 'onerow');
axh = axes('Position', [0 0.25 .3 .7]);  % Position vector: [left bottom width height]
gencue_disp = montage(gencue_disp, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);
gencue_disp = removeblobs(gencue_disp)

% positive
gencue_disp = addblobs(gencue_disp, clpos_data_gencue_A{1}, 'color', [.8 .2 0]);  % path a 
gencue_disp = addblobs(gencue_disp, clpos_data_gencue_B{1}, 'outline', [.4 .4 .4]);  % path b 
gencue_disp = addblobs(gencue_disp, clpos_data_gencue_AB{1}, 'color', [0 1 0]);  % path ab effect

% % negative
% gencue_disp = addblobs(gencue_disp, clneg_data_gencue_A{1}, 'color', [0 .1 .6]);  
% gencue_disp = addblobs(gencue_disp, clneg_dat_gencue_B{1}, 'color', [0 .4 .8]);  
% gencue_disp = addblobs(gencue_disp, clneg_data_gencue_AB{1}, 'color', [0 .8 .8]); 


%% (FDR corrected) mediation results during PAIN

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/mediation_38/SimHC_WBgenpain_painrate_36

SETUP = mediation_brain_corrected_threshold('fdr', 'mask', 'mask.img');  % to get FDR threshold across all mediation result images  

[clpos_genpain_A_fdr, clneg_genpain_A_fdr, clpos_data_genpain_A_fdr, clneg_data_genpain_A_fdr] = mediation_brain_results('a', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask)
[clpos_genpain_A, clneg_genpain_A, clpos_data_genpain_A, clneg_data_genpain_A] = mediation_brain_results('a', 'thresh', [SETUP.fdr_p_thresh .01], 'size', [3 1], 'mask', mask, 'prune')
[clpos_genpain_B_fdr, clneg_genpain_B_fdr, clpos_data_genpain_B_fdr, clneg_dat_genpain_B_fdr] = mediation_brain_results('b', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask)
[clpos_genpain_B, clneg_genpain_B, clpos_data_genpain_B, clneg_dat_genpain_B] = mediation_brain_results('b', 'thresh', [SETUP.fdr_p_thresh .01], 'size', [3 1], 'mask', mask, 'prune')
[clpos_genpain_AB_fdr, clneg_genpain_AB_fdr, clpos_data_genpain_AB_fdr, clneg_data_genpain_AB_fdr] = mediation_brain_results('ab', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask)
[clpos_genpain_AB, clneg_genpain_AB, clpos_data_genpain_AB, clneg_data_genpain_AB] = mediation_brain_results('ab', 'thresh', [SETUP.fdr_p_thresh .01], 'size', [3 1], 'mask', mask, 'prune')

xmy_genpain = fmri_data('X-M-Y_effect.img');


% close all

%% make plots for PAIN mediators

genpain_disp = fmridisplay;
% genpain_disp = montage(genpain_disp, 'axial', 'wh_slice', [0 0 -66; 0 0 -10; 0 0 10; 0 0 30; 0 0 40; 0 0 50], 'onerow');
genpain_disp = montage(genpain_disp, 'axial', 'wh_slice', [0 0 -66; 0 0 -10; 0 0 0; 0 0 10; 0 0 20; 0 0 30; 0 0 40; 0 0 50; 0 0 60], 'onerow');
axh = axes('Position', [0 0.25 .3 .7]);  % Position vector: [left bottom width height]
genpain_disp = montage(genpain_disp, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh);
genpain_disp = removeblobs(genpain_disp)

% positive
genpain_disp = addblobs(genpain_disp, clpos_data_genpain_A{1}, 'color', [.8 .2 0]);  % path a 
genpain_disp = addblobs(genpain_disp, clpos_data_genpain_B{1}, 'outline', [.6 .4 0]);  % path b 
genpain_disp = addblobs(genpain_disp, clpos_data_genpain_AB{1}, 'color', [.9 .9 0]);  % path ab effect


