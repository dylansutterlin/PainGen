%% save clusters to images (to use them later or plot)

clear all
close all
clc

%% mediation results

cd C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36

SETUP = mediation_brain_corrected_threshold('fdr', 'mask', 'mask.img');  % to get FDR threshold across all mediation result images  

mask = which('gray_matter_mask.img');

[clpos_gen_A_fdr, clneg_gen_A_fdr, clpos_data_gen_A_fdr, clneg_data_gen_A_fdr] = mediation_brain_results('a', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask, 'noplots');
[clpos_gen_A, clneg_gen_A, clpos_data_gen_A, clneg_data_gen_A] = mediation_brain_results('a', 'thresh', [SETUP.fdr_p_thresh .01], 'size', [3 1], 'mask', mask, 'prune', 'noplots');
[clpos_gen_B_fdr, clneg_gen_B_fdr, clpos_data_gen_B_fdr, clneg_data_gen_B_fdr] = mediation_brain_results('b', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask, 'noplots');
[clpos_gen_B, clneg_gen_B, clpos_data_gen_B, clneg_data_gen_B] = mediation_brain_results('b', 'thresh', [SETUP.fdr_p_thresh .01], 'size', [3 1], 'mask', mask, 'prune', 'noplots');
[clpos_gen_AB_fdr, clneg_gen_AB_fdr, clpos_data_gen_AB_fdr, clneg_data_gen_AB_fdr] = mediation_brain_results('ab', 'thresh', SETUP.fdr_p_thresh, 'size', 3, 'mask', mask, 'noplots');
[clpos_gen_AB, clneg_gen_AB, clpos_data_gen_AB, clneg_data_gen_AB] = mediation_brain_results('ab', 'thresh', [SETUP.fdr_p_thresh .01], 'size', [3 1], 'mask', mask, 'prune', 'noplots');

%% hippocampus AB cluster mask

% cluster #15
hippAB = cluster2region(clpos_data_gen_AB_fdr{1, 1}(15));
orthviews(hippAB)
hippo_mask = region2imagevec(hippAB);
hippo_mask.fullpath = 'C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results\hippocampus_AB_mask.img';
write(hippo_mask)

A_hippomask = apply_mask(fmri_data('X-M_indiv_effect.img'), hippo_mask.fullpath);
A_hippo_mean = mean(A_hippomask.dat);

AB_hippomask = apply_mask(fmri_data('X-M-Y_indiv_effect.img'), hippo_mask.fullpath);
AB_hippo_mean = mean(AB_hippomask.dat);

figure; subplot(1,2,1); scatter(cue_exp_effect_N36, AB_hippo_mean); refline
subplot(1,2,2); scatter(cue_exp_effect_N36, A_hippo_mean); refline

%% write CSCat results

cd C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results

genAr = cluster2region(clpos_data_gen_A_fdr{1});
Gen_A_fdr = region2imagevec(genAr);
Gen_A_fdr.fullpath = 'C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results\Gen_Apos_fdr.img';
write(Gen_A_fdr)

genABr = cluster2region(clpos_data_gen_AB_fdr{1});
Gen_AB_fdr = region2imagevec(genABr);
Gen_AB_fdr.fullpath = 'C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results\Gen_ABpos_fdr.img';
write(Gen_AB_fdr)

genBr = cluster2region(clpos_data_gen_B_fdr{1});
Gen_B_fdr = region2imagevec(genBr);
Gen_B_fdr.fullpath = 'C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results\Gen_Bpos_fdr.img';
write(Gen_B_fdr)

%% for 01 (pruned from FDR)

genAr = cluster2region(clpos_data_gen_A{1});
Gen_A_01 = region2imagevec(genAr);
Gen_A_01.fullpath = 'C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results\Gen_Apos_01.img';
write(Gen_A_01)

genABr = cluster2region(clpos_data_gen_AB{1});
Gen_AB_01 = region2imagevec(genABr);
Gen_AB_01.fullpath = 'C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results\Gen_ABpos_01.img';
write(Gen_AB_01)

genBr = cluster2region(clpos_data_gen_B{1});
Gen_B_01 = region2imagevec(genBr);
Gen_B_01.fullpath = 'C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Imaging\Blanca_2017\180710_CSCat_WBgenpain_painrate_36\Results\Gen_Bpos_01.img';
write(Gen_B_01)



