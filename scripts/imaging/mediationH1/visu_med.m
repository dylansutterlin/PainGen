%% 

resultDir = '/home/dsutterlin/projects/genPain/results/imaging/mediation/'

%% Cue-M-Pain model
mask = which('gray_matter_mask.img');
%SETUP = mediation_brain_corrected_threshold('fdr');
SETUP = mediation_brain_corrected_threshold('fdr', 'mask', 'mask.img');  % to get FDR threshold across all mediation result images  

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('\n\n%s\n%s\n%s\n%s\n%s\n\n', dashes, dashes, str, dashes, dashes);

printhdr('Path a')
mediation_brain_results('a', 'thresh', ...
    SETUP.fdr_p_thresh, 'size', 1, ...
    'slices', 'tables', 'names', 'save');
mediation_brain_results('a', 'thresh', [0.002 .05], 'size', [1 5], 'mask', mask, 'slices', 'tables', 'names', 'save')

printhdr('Path b:')
mediation_brain_results('b', 'thresh', ...
    [SETUP.fdr_p_thresh .05], 'size', [1 5], ...
    'slices', 'tables', 'names', 'save');

printhdr('Path a*b:')
[clpos, clneg, clpos_data, clneg_data, clpos_data2, clneg_data2] = mediation_brain_results('ab', 'thresh', ...
    [SETUP.fdr_p_thresh .05], 'size', [1 5], ...
    'slices', 'tables', 'names', 'save');

cluster_table(clpos_data, 1, 0, 'num_sig_voxels')
 

[clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all', 'thresh', [Inf .005 .01], 'size', [1 3 10], 'fdrthresh', .05, 'prune', 'conj');
[clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all', 'thresh', [.005 .05 .05], 'size', [1 5 10], 'fdrthresh', .05, 'prune');
iv = InteractiveViewer(args{:}, 'UseExistingGraphicsWindow', 1, 'LoadDataOnDemand', 1, 'IVObserver', mediationIVObserver('X', X, 'Y', Y));

[clpos_gencue_A_001, clneg_gencue_A_001, clpos_data_gencue_A_001, clneg_data_gencue_A_001] = mediation_brain_results('a', 'thresh', .001, 'size', 5, 'mask', mask)
[clpos_gencue_A, clneg_gencue_A, clpos_data_gencue_A, clneg_data_gencue_A] = mediation_brain_results('a', 'thresh', [.001 .01], 'size', [5 20], 'mask', mask)
[clpos_gencue_B_001, clneg_gencue_B_001, clpos_data_gencue_B_001, clneg_data_gencue_B_001] = mediation_brain_results('b', 'thresh', .001, 'size', 5, 'mask', mask)
[clpos_gencue_B, clneg_gencue_B, clpos_data_gencue_B, clneg_data_gencue_B] = mediation_brain_results('b', 'thresh', [.001 .01], 'size', [5 20], 'mask', mask)
[clpos_gencue_AB_001, clneg_gencue_AB_001, clpos_data_gencue_AB_001, clneg_data_gencue_AB_001] = mediation_brain_results('ab', 'thresh', .001, 'size', 5, 'mask', mask)
[clpos_gencue_AB, clneg_gencue_AB, clpos_data_gencue_AB, clneg_data_gencue_AB] = mediation_brain_results('ab', 'thresh', [.001 .01], 'size', [5 20], 'mask', mask)

