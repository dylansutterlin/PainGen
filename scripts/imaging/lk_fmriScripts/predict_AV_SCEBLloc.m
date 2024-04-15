%% try to develop classifier separating animals from vehicles

clear all
close all
clc

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/model4_loc

overlay = '/Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/mean_SCEBL_anatomical_sharp_noskull.img';


%% contrasts 

Acons = filenames('SCEBL*/con_0013.img');
Vcons = filenames('SCEBL*/con_0014.img');
Hcons = filenames('SCEBL*/con_0015.img');
Tcons = filenames('SCEBL*/con_0016.img');

%% put into fmri_data object
AV_allcons = fmri_data([Acons; Vcons], which('gray_matter_mask.img'));
plot(AV_allcons)
AV_allcons.Y = [ones(37,1); -1*ones(37,1)];

AVHT_allcons = fmri_data([Acons; Vcons; Tcons; Hcons], (Acons{1}));
plot(AVHT_allcons)

% sort by living(natural) vs artificial
AVHT_allcons.Y = [ones(37,1); -1*ones(37,1); ones(37,1); -1*ones(37,1)];

% sort by moving vs static
AVHT_allcons.Y = [ones(37,1); ones(37,1); -1*ones(37,1); -1*ones(37,1)];


%% single exemplar contrasts for each of the four categories

A1 = fmri_data(filenames('SCEBL*/con_0001.img'));
A2 = fmri_data(filenames('SCEBL*/con_0002.img'));
A3 = fmri_data(filenames('SCEBL*/con_0003.img'));
V1 = fmri_data(filenames('SCEBL*/con_0004.img'));
V2 = fmri_data(filenames('SCEBL*/con_0005.img'));
V3 = fmri_data(filenames('SCEBL*/con_0006.img'));
H1 = fmri_data(filenames('SCEBL*/con_0007.img'));
H2 = fmri_data(filenames('SCEBL*/con_0008.img'));
H3 = fmri_data(filenames('SCEBL*/con_0009.img'));
T1 = fmri_data(filenames('SCEBL*/con_0010.img'));
T2 = fmri_data(filenames('SCEBL*/con_0011.img'));
T3 = fmri_data(filenames('SCEBL*/con_0012.img'));

A123_VTH = fmri_data([filenames('SCEBL*/con_0001.img'); filenames('SCEBL*/con_0002.img'); filenames('SCEBL*/con_0003.img'); filenames('SCEBL*/con_0014.img'); filenames('SCEBL*/con_0016.img'); filenames('SCEBL*/con_0015.img')]);
V123_ATH = fmri_data([filenames('SCEBL*/con_0004.img'); filenames('SCEBL*/con_0005.img'); filenames('SCEBL*/con_0006.img'); filenames('SCEBL*/con_0013.img'); filenames('SCEBL*/con_0016.img'); filenames('SCEBL*/con_0015.img')]);
% T123 = fmri_data([filenames('SCEBL*/con_0010.img'); filenames('SCEBL*/con_0011.img'); filenames('SCEBL*/con_00012.img')]);
% H123 = fmri_data([filenames('SCEBL*/con_0007.img'); filenames('SCEBL*/con_0008.img'); filenames('SCEBL*/con_0009.img')]);

A123_V123 = fmri_data([filenames('SCEBL*/con_0001.img'); filenames('SCEBL*/con_0002.img'); filenames('SCEBL*/con_0003.img'); filenames('SCEBL*/con_0004.img'); filenames('SCEBL*/con_0005.img'); filenames('SCEBL*/con_0006.img')]);

A123_VTH.Y = [ones(37,1); ones(37,1); ones(37,1); -1*ones(37,1); -1*ones(37,1); -1*ones(37,1)];
V123_ATH.Y = [ones(37,1); ones(37,1); ones(37,1); -1*ones(37,1); -1*ones(37,1); -1*ones(37,1)];
A123_V123.Y = [ones(37,1); ones(37,1); ones(37,1); -1*ones(37,1); -1*ones(37,1); -1*ones(37,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% predict

addpath(genpath('/Users/Leonie/Documents/Tools/spider'))

[cverr, stats, optout] = predict(AV_allcons, 'algorithm_name', 'cv_svm', 'nfolds', [1:37,1:37]);

% [cverr, stats, optout] = predict(A123_VTH, 'algorithm_name', 'cv_svm', 'nfolds', 'loso', 'error_type', 'mse');
% 
% [cverr, stats, optout] = predict(V123_ATH, 'algorithm_name', 'cv_svm', 'nfolds', 'loso', 'error_type', 'mse');
% 
% [cverr, stats, optout] = predict(A123_V123, 'algorithm_name', 'cv_svm', 'nfolds', 'loso', 'error_type', 'mse');
% 

%% threshold weight object arbitrarily and plot it

AVmap = fmri_data(stats.weight_obj);
orthviews(AVmap);

% save map
cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/model4_loc/AV_predict_map
AVmap.fullpath = fullfile(pwd, 'AVmap.img');
write(AVmap)

% close all

%% z-score and plot map
AVmapz = rescale(AVmap, 'zscoreimages')
% plot(AVmapz)
% orthviews(AVmapz)

trAVmapz = threshold(AVmapz, [-2 2], 'raw-outside');
AVregs = region(trAVmapz)
% orthviews(trAVmapz)
% surface(trgmapz)

AVmap_disp = fmridisplay(AVmap) %, 'overlay', overlay)
AVmap_disp = montage(AVmap_disp, 'saggital', 'slice_range', [-50 50],  'spacing', 10);
AVmap_disp = addblobs(AVmap_disp, AVregs)
AVmap_disp = montage(AVmap_disp, 'coronal', 'slice_range', [-50 50],  'spacing', 10);
AVmap_disp = addblobs(AVmap_disp, AVregs)
AVmap_disp = montage(AVmap_disp, 'transversal', 'slice_range', [-40 60],  'spacing', 5);
AVmap_disp = addblobs(AVmap_disp, AVregs)


%% ROC Plots 

create_figure('AnimalVehicle_Twochoice_ROC')
ROC_twochoice = roc_plot(stats.dist_from_hyperplane_xval, stats.Y == 1, 'twochoice', 'color', 'r');
set(gca, 'OuterPosition', [0 0 1 1]);


%% apply classifier to all of them

pexpvals_A1 = apply_mask(A1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_A2 = apply_mask(A2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_A3 = apply_mask(A3, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_V1 = apply_mask(V1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_V2 = apply_mask(V2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_V3 = apply_mask(V3, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_T1 = apply_mask(T1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_T2 = apply_mask(T2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_T3 = apply_mask(T3, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_H1 = apply_mask(H1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_H2 = apply_mask(H2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_H3 = apply_mask(H3, AVmap, 'pattern_expression', 'ignore_mising'); 

% concatenate them

pexpvals_all = [pexpvals_A1 pexpvals_A2 pexpvals_A3 ...
                pexpvals_V1 pexpvals_V2 pexpvals_V3 ...
                pexpvals_T1 pexpvals_T2 pexpvals_T3 ...
                pexpvals_H1 pexpvals_H2 pexpvals_H3];

%% Figure

create_figure('Categories_AVPatternexp', 1, 2);
barplot_colored(pexpvals_all, 'Pattern Expression', [], 'y',(-1:1)+.5,'nofig', 'noind', 'within');
set(gca, 'XTick',[1:12], 'XTickLabel',{'', 'Animals', '', '', 'Vehicles', '', '', 'Trees', '', '', 'Houses', '', });
cm = [1 .5 0; .5 0 1; 0 .5 .2; 0 .1 .8];
colormap(cm)
xlabel('Session Conditions');
ylabel('PExp Values');

subplot(1, 2, 2)
plot(pexpvals_all(dat.Y > 0), 'o','color',[1 .5 0], 'LineWidth', 3)
plot(pexpvals_all(dat.Y < 0), 'o','color',[.5 0 1], 'LineWidth', 3)
ylabel('PExp Values');
xlabel('Subject Number')



%% use on main task cue representation (and see whether this predicts learning?)

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/model4_34

AV_h = [2, 1, 2, 1, 2, ...
        0, 2, 1, 2, 0, ...
        1, 2, 1, 2, 1, ...
        2, 1, 2, ...
        2, 1, 2, 1, 2, ...
        1, 2, 1, 2, 1, ...
        1, 2, 1, 2, 1, ...
        2, 1, 1, 1];  % missing 236

CLSL = fmri_data(filenames('SCEBL*/con_0001.img'));
CLSH = fmri_data(filenames('SCEBL*/con_0002.img'));
CHSL = fmri_data(filenames('SCEBL*/con_0003.img'));
CHSH = fmri_data(filenames('SCEBL*/con_0004.img'));

Pain49CLSL = fmri_data(filenames('SCEBL*/con_0008.img'));
Pain49CLSH = fmri_data(filenames('SCEBL*/con_0009.img'));
Pain49CHSL = fmri_data(filenames('SCEBL*/con_0010.img'));
Pain49CHSH = fmri_data(filenames('SCEBL*/con_0011.img'));

CH_CL = fmri_data(filenames('SCEBL*/con_0020.img'));

pexpvals_CLSL = apply_mask(CLSL, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 
pexpvals_CLSH = apply_mask(CLSH, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 
pexpvals_CHSL = apply_mask(CHSL, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 
pexpvals_CHSH = apply_mask(CHSH, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 

pexpvals_Pain49CLSL = apply_mask(Pain49CLSL, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 
pexpvals_Pain49CLSH = apply_mask(Pain49CLSH, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 
pexpvals_Pain49CHSL = apply_mask(Pain49CHSL, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 
pexpvals_Pain49CHSH = apply_mask(Pain49CHSH, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 


pexpvals_CH_CL = apply_mask(CH_CL, stats.weight_obj, 'pattern_expression', 'ignore_mising'); 

pexpvals_learning = [pexpvals_CLSL pexpvals_CLSH pexpvals_CHSL pexpvals_CHSH];

% separated by counterbalancing conditions (A or V high cue)
pexpvals_learning_A_CH = pexpvals_learning(AV_h == 1, :);
pexpvals_learning_V_CH = pexpvals_learning(AV_h == 2, :);

create_figure('Categories_AVPatternexp LEARNING by CS low vs high');
barplot_colored(pexpvals_learning, 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within');
set(gca, 'XTick',[1:4], 'XTickLabel',{'CLSL', 'CLSH', 'CHSL', 'CHSH'});
cm = [1 .5 0; .5 0 1; 0 .5 .2; 0 .1 .8];
colormap(cm)
xlabel('Learning Conditions');
ylabel('PExp Values');

create_figure('Categories_AVPatternexp LEARNING by CS low vs high, separated for A and V CH');
subplot(1,2,1)
barplot_colored(pexpvals_learning_A_CH, 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within'); hold on
plot(1:4, pexpvals_learning_A_CH, 'ok')
set(gca, 'XTick', 1:4, 'XTickLabel', {'CLSL', 'CLSH', 'CHSL', 'CHSH'});
xlabel('Learning Conditions CShigh = Animal');
ylabel('PExp Values');
subplot(1,2,2)
barplot_colored(pexpvals_learning_V_CH, 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within'); hold on
plot(1:4, pexpvals_learning_V_CH, 'ok')
set(gca, 'XTick', 1:4, 'XTickLabel', {'CLSL', 'CLSH', 'CHSL', 'CHSH'});
colormap(cm)
xlabel('Learning Conditions CShigh = Vehicle');
ylabel('PExp Values');
    
create_figure('Categories_AVPatternexp LEARNING CH > CL, separated for A and V CH');
subplot(1,2,1)
barplot_colored(pexpvals_CH_CL(AV_h == 1, :), 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within'); hold on
plot(1, pexpvals_CH_CL(AV_h == 1, :), 'ok')
xlabel('Learning Conditions CShigh = Animal');
ylabel('PExp Values');
subplot(1,2,2)
barplot_colored(pexpvals_CH_CL(AV_h == 2, :), 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within'); hold on
plot(1, pexpvals_CH_CL(AV_h == 2, :), 'ok')
xlabel('Learning Conditions CShigh = Vehicle');
ylabel('PExp Values');

pexpvals_CH_CL_rev = [pexpvals_CH_CL(AV_h == 1); -1* pexpvals_CH_CL(AV_h == 2)];

[tmain pmain] = ttest([nanmean(pexpvals_learning_A_CH(:,1:2), 2); nanmean(pexpvals_learning_V_CH(:,3:4), 2)], [nanmean(pexpvals_learning_A_CH(:,3:4), 2); nanmean(pexpvals_learning_V_CH(:,1:2), 2)]);


%% during pain

pexpvals_learning_Pain49 = [pexpvals_Pain49CLSL pexpvals_Pain49CLSH pexpvals_Pain49CHSL pexpvals_Pain49CHSH];

create_figure('Categories_AVPatternexp LEARNING during PAIN (49), separated for A and V CH');
subplot(1,2,1)
barplot_colored(pexpvals_learning_Pain49(AV_h == 1, :), 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within'); hold on
plot(1:4, pexpvals_learning_Pain49(AV_h == 1, :), 'ok')
set(gca, 'XTick', 1:4, 'XTickLabel', {'CLSL', 'CLSH', 'CHSL', 'CHSH'});
xlabel('Learning Conditions CShigh = Animal');
ylabel('PExp Values');
subplot(1,2,2)
barplot_colored(pexpvals_learning_Pain49(AV_h == 2, :), 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within'); hold on
plot(1:4, pexpvals_learning_Pain49(AV_h == 2, :), 'ok')
set(gca, 'XTick', 1:4, 'XTickLabel', {'CLSL', 'CLSH', 'CHSL', 'CHSH'});
colormap(cm)
xlabel('Learning Conditions CShigh = Vehicle');
ylabel('PExp Values');

create_figure('Categories_AVPatternexp PAIN');
barplot_colored(pexpvals_learning_Pain49, 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within');
set(gca, 'XTick',[1:4], 'XTickLabel',{'CLSL', 'CLSH', 'CHSL', 'CHSH'});
cm = [1 .5 0; .5 0 1; 0 .5 .2; 0 .1 .8];
colormap(cm)
xlabel('Learning Conditions');
ylabel('PExp Values');




%% test expression in generalization task

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/model2_gen

GCA1 = fmri_data(filenames('SCEBL*/con_0001.img'));
GCA2 = fmri_data(filenames('SCEBL*/con_0002.img'));
GCA3 = fmri_data(filenames('SCEBL*/con_0003.img'));
GPA1 = fmri_data(filenames('SCEBL*/con_0004.img'));
GPA2 = fmri_data(filenames('SCEBL*/con_0005.img'));
GPA3 = fmri_data(filenames('SCEBL*/con_0006.img'));
GWA1 = fmri_data(filenames('SCEBL*/con_0007.img'));
GWA2 = fmri_data(filenames('SCEBL*/con_0008.img'));
GWA3 = fmri_data(filenames('SCEBL*/con_0009.img'));

GCV1 = fmri_data(filenames('SCEBL*/con_0010.img'));
GCV2 = fmri_data(filenames('SCEBL*/con_0011.img'));
GCV3 = fmri_data(filenames('SCEBL*/con_0012.img'));
GPV1 = fmri_data(filenames('SCEBL*/con_0013.img'));
GPV2 = fmri_data(filenames('SCEBL*/con_0014.img'));
GPV3 = fmri_data(filenames('SCEBL*/con_0015.img'));
GWV1 = fmri_data(filenames('SCEBL*/con_0016.img'));
GWV2 = fmri_data(filenames('SCEBL*/con_0017.img'));
GWV3 = fmri_data(filenames('SCEBL*/con_0018.img'));


pexpvals_GCA1 = apply_mask(GCA1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GCA2 = apply_mask(GCA2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GCA3 = apply_mask(GCA3, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GPA1 = apply_mask(GPA1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GPA2 = apply_mask(GPA2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GPA3 = apply_mask(GPA3, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GWA1 = apply_mask(GWA1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GWA2 = apply_mask(GWA2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GWA3 = apply_mask(GWA3, AVmap, 'pattern_expression', 'ignore_mising'); 

pexpvals_GCV1 = apply_mask(GCV1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GCV2 = apply_mask(GCV2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GCV3 = apply_mask(GCV3, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GPV1 = apply_mask(GPV1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GPV2 = apply_mask(GPV2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GPV3 = apply_mask(GPV3, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GWV1 = apply_mask(GWV1, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GWV2 = apply_mask(GWV2, AVmap, 'pattern_expression', 'ignore_mising'); 
pexpvals_GWV3 = apply_mask(GWV3, AVmap, 'pattern_expression', 'ignore_mising'); 

pexpvals_gen = [pexpvals_GCA1 pexpvals_GCA2 pexpvals_GCA3 ...
                pexpvals_GPA1 pexpvals_GPA2 pexpvals_GPA3 ...
                pexpvals_GWA1 pexpvals_GWA2 pexpvals_GWA3 ...
                pexpvals_GCV1 pexpvals_GCV2 pexpvals_GCV3 ...
                pexpvals_GPV1 pexpvals_GPV2 pexpvals_GPV3 ...
                pexpvals_GWV1 pexpvals_GWV2 pexpvals_GWV3];
            
%% plot            

cmgen = [0 .9 .3; 0 .8 .2; 0 .7 .1; 0 .8 .6; 0 .7 .5; 0 .6 .4; 0 .7 .9; 0 .6 .8; 0 .5 .7; ...
         .9 0 .3; .8 0 .2; .7 0 .1; .8 0 .6; .7 0 .5; .6 0 .4; .7 0 .9; .6 0 .8; .5 0 .7];

create_figure('Categories_AVPatternexp GENERAL');
barplot_colored(pexpvals_gen, 'Pattern Expression', [], 'y',(-1:1)+.5,'nofig', 'noind', 'within');
set(gca, 'XTick',[1:18], 'XTickLabel',{'', 'ANI cartoons', '', '', 'ANI photos', '', '', 'ANI words', '', '', 'VEH cartoons', '', '', 'VEH photos', '','', 'VEH words', ''});
cm = [1 .5 0; .5 0 1; 0 .5 .2; 0 .1 .8];
colormap(cmgen)
xlabel('Session Conditions');
ylabel('PExp Values');






