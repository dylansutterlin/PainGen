% Thresholded results displays
% April 2018 LK, based on many previous versions of this script

clear 
close all
clc

a_cols = {[.9 .4 0]; [1 .6 0]};
ab_cols = {[.7 .2 .3]; [1 .4 .55]};

% % other color options? more mint green and peach red
% so_cols = {[.7 .2 .3]; [1 .4 .55]};
% cs_cols = {[.3 .6 .2]; [.5 .9 .45]};


%% Have both in the same plots, by path

mask = which('gray_matter_mask.img');

% SETUP = mediation_brain_corrected_threshold('fdr', 'mask', 'mask.img');  % to get FDR threshold across all mediation result images  
% [clpos_A_fdr, clneg_A_fdr, clpos_data_A_fdr, clneg_data_A_fdr] = mediation_brain_results('a', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask, 'noplots');
% [clpos_A, clneg_A, clpos_data_A, clneg_data_A] = mediation_brain_results('a', 'thresh', [SETUP.fdr_p_thresh .01 .05], 'size', [3 1 1], 'mask', mask, 'prune', 'noplots');
% [clpos_B_fdr, clneg_B_fdr, clpos_data_B_fdr, clneg_data_B_fdr] = mediation_brain_results('b', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask, 'noplots');
% [clpos_B, clneg_B, clpos_data_B, clneg_data_B] = mediation_brain_results('b', 'thresh', [SETUP.fdr_p_thresh .01 .05], 'size', [3 1 1], 'mask', mask, 'prune', 'noplots');
% [clpos_AB_fdr, clneg_AB_fdr, clpos_data_AB_fdr, clneg_data_AB_fdr] = mediation_brain_results('ab', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask, 'noplots');
% [clpos_AB, clneg_AB, clpos_data_AB, clneg_data_AB] = mediation_brain_results('ab', 'thresh', [SETUP.fdr_p_thresh .01 .05], 'size', [3 1 1], 'mask', mask, 'prune', 'noplots');


%%

cd C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results

Gen_A_pos01 = fmri_data('Gen_Apos_01.img');
Gen_A_posfdr = fmri_data('Gen_Apos_fdr.img');
Gen_B_pos01 = fmri_data('Gen_Bpos_01.img');
Gen_B_posfdr = fmri_data('Gen_Bpos_fdr.img');
Gen_AB_pos01 = fmri_data('Gen_ABpos_01.img');
Gen_AB_posfdr = fmri_data('Gen_ABpos_fdr.img');
% CS_conj = fmri_data('intersection_A_AB_05pruned.img');


%% display

disp_A = fmridisplay;
disp_A = montage(disp_A, 'axial', 'wh_slice', [0 0 -22; 0 0 -24; 0 0 -12], 'onerow');
% axh = axes('Position', [0 0.25 .3 .7]);  % Position vector: [left bottom width height]
% disp_A = montage(disp_A, 'saggital', 'wh_slice', [-4 0 0], 'existing_axes', axh); title('Path a')
disp_A = montage(disp_A, 'saggital', 'wh_slice', [0 0 0; 0 0 0; -26 0 0; -24 0 0], 'onerow'); % hippocampus and some amygdala
disp_A = montage(disp_A, 'coronal', 'wh_slice', [0 0 0; 0 0 0; 0 0 0]); % amydala

disp_A = removeblobs(disp_A);
disp_A = addblobs(disp_A, region(Gen_A_pos01), 'color', a_cols{1});  % path a
disp_A = addblobs(disp_A, region(Gen_A_posfdr), 'color', a_cols{2});  % path a

%% for B path

disp_B = fmridisplay;
disp_B = montage(disp_B, 'axial', 'wh_slice', [0 0 -66; 0 0 0; 0 0 20; 0 0 40], 'onerow');
axh = axes('Position', [0 0.25 .3 .7]);  % Position vector: [left bottom width height]
disp_B = montage(disp_B, 'saggital', 'wh_slice', [-6 0 0], 'existing_axes', axh); title('Path b')
disp_B = removeblobs(disp_B);

disp_B = addblobs(disp_B, region(Gen_B_pos01), 'color', [.8 .5 .08]);  % path b 
disp_B = addblobs(disp_B, region(Gen_B_posfdr), 'color', [1 .8 .2]);  % path b 

% size
set(gcf, 'Position', [100 100 1200 400])


%% for AB path

disp_AB = fmridisplay;
disp_AB = montage(disp_AB, 'axial', 'wh_slice', [0 0 -28; 0 0 -26; 0 0 -24; 0 0 -16; 0 0 -10; 0 0 -8], 'onerow');
disp_AB = montage(disp_AB, 'axial', 'wh_slice', [0 0 -26; 0 0 -10; 0 0 -6], 'onerow');

disp_AB = montage(disp_AB, 'saggital', 'wh_slice', [-32 0 0; -28 0 0; 2 0 0], 'onerow'); % -31 ist best here for hippocampus and amygdala

% disp_AB = montage(disp_AB, 'coronal', 'wh_slice', [0 6 0; 0 2 0; 0 0 0; 0 -4 0], 'onerow'); % -18 is good, -10 (bilat) and -6 (right) for ant hippoc

disp_AB = removeblobs(disp_AB);

disp_AB = addblobs(disp_AB, region(Gen_AB_pos01), 'color', ab_cols{1});  % path a 
disp_AB = addblobs(disp_AB, region(Gen_AB_posfdr), 'color', ab_cols{2});  % path a 

