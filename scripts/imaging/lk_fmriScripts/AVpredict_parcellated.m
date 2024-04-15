% by LK June 2015

clear all
close all
clc

addpath(genpath('/Users/Leonie/Documents/Tools/spider'))

%% load data

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/model4_loc
Acons = filenames('SCEBL*/con_0013.img');
Vcons = filenames('SCEBL*/con_0014.img');
Hcons = filenames('SCEBL*/con_0015.img');
Tcons = filenames('SCEBL*/con_0016.img');

AV_allcons = fmri_data([Acons; Vcons], which('gray_matter_mask.img'));
% plot(AV_allcons)
AV_allcons.Y = [ones(37,1); -1*ones(37,1)];


%% get and plot parcellation data (creating region object and masks)

% % Using 4297 parcels from Braimap (?)
% parcels = fmri_data(which('brainmap_july2010_4297_parcels.img'));

% Using Cameron Craddock's ADHD sample parcels (190) based on functional parcellation
parcels = fmri_data('/Users/Leonie/Documents/Tools/ADHD200_parcellations/ADHD200_parcellate_200.nii')
% plot(parcels)
parcel_regs = region(parcels, 'unique_mask_values');
orthviews(parcel_regs)

for r = [1:34 36:44 46:numel(parcel_regs)]
    clear mask maskedAV
    mask = region2imagevec(parcel_regs(r));         % create mask from region object
    maskedAV = apply_mask(AV_allcons, mask);       % apply mask to fmri_data object
    % orthviews(maskedAV)
        
    [cverr{r}, stats{r}, optout{r}] = predict(maskedAV, 'algorithm_name', 'cv_svm', 'nfolds', [1:37,1:37]);
    ROC_twochoice{r} = roc_plot(stats{r}.dist_from_hyperplane_xval, stats{r}.Y == 1, 'twochoice', 'color', 'r'); close all
    
    if ROC_twochoice{r}.accuracy_p < 0.05
        sign05(r) = 1;
    else sign05(r) = 0;
    end
       
end


%% now show all significant regions

figure; hist(cell2mat(cverr))

for r = [1:34 36:44 46:numel(parcel_regs)]
   
    if ROC_twochoice{r}.accuracy_p < 0.01
        sign01(r) = 1;
        if ROC_twochoice{r}.accuracy_p < 0.001
            sign001(r) = 1;
        else sign001(r) = 0;
        end
    else sign01(r) = 0;
    end
       
end

signparcels = parcel_regs(find(sign01));
orthviews(signparcels)


%% use those to create new maska and bootstrapped patterns?



%% use on main task cue representation (and see whether this predicts learning)

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/model4_34

AV_h = [2, 1, 2, 1, 2,...
        1, 2, 1, 2, 1,  ...
        1, 2, 1, 2, 1, ...
        2, 1, 2, ...
        2, 1, 2, 1, 2,...
        1, 2, 1, 2, 1,...
        1, 2, 1, 2, 1,...
        2, 1, 1, 1];  % missing 236

CLSL = fmri_data(filenames('SCEBL*/con_0001.img'));
CLSH = fmri_data(filenames('SCEBL*/con_0002.img'));
CHSL = fmri_data(filenames('SCEBL*/con_0003.img'));
CHSH = fmri_data(filenames('SCEBL*/con_0004.img'));

Pain49CLSL = fmri_data(filenames('SCEBL*/con_0008.img'));
Pain49CLSH = fmri_data(filenames('SCEBL*/con_0009.img'));
Pain49CHSL = fmri_data(filenames('SCEBL*/con_0010.img'));
Pain49CHSH = fmri_data(filenames('SCEBL*/con_0011.img'));


for sigr = 1:numel(find(sign01))
    
    pexpvals_CLSL = apply_mask(CLSL, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_CLSH = apply_mask(CLSH, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_CHSL = apply_mask(CHSL, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_CHSH = apply_mask(CHSH, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 

    pexpvals_learning{sigr} = [pexpvals_CLSL pexpvals_CLSH pexpvals_CHSL pexpvals_CHSH];
    
    % separated by counterbalancing conditions (A or V high cue)
    pexpvals_learning_A_CH{sigr} = pexpvals_learning{sigr}(AV_h == 1, :);
    pexpvals_learning_V_CH{sigr} = pexpvals_learning{sigr}(AV_h == 2, :);
    
    clear pexpvals_C*
    
    create_figure('Categories_AVPatternexp LEARNING by CS low vs high');
    barplot_colored(pexpvals_learning{sigr}, 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within');
    set(gca, 'XTick',[1:4], 'XTickLabel',{'CLSL', 'CLSH', 'CHSL', 'CHSH'});
    cm = [1 .5 0; .5 0 1; 0 .5 .2; 0 .1 .8];
    colormap(cm)
    xlabel('Learning Conditions');
    ylabel('PExp Values');
    
    create_figure('Categories_AVPatternexp LEARNING by CS low vs high, separated for A and V CH');
    subplot(1,2,1)
    barplot_colored(pexpvals_learning_A_CH{sigr}, 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within');
    set(gca, 'XTick', 1:4, 'XTickLabel', {'CLSL', 'CLSH', 'CHSL', 'CHSH'});
    xlabel('Learning Conditions CShigh = Animal');
    ylabel('PExp Values');
    subplot(1,2,2)
    barplot_colored(pexpvals_learning_V_CH{sigr}, 'Pattern Expression', [], 'y', (-1:1)+.5, 'nofig', 'noind', 'within');
    set(gca, 'XTick', 1:4, 'XTickLabel', {'CLSL', 'CLSH', 'CHSL', 'CHSH'});
    colormap(cm)
    xlabel('Learning Conditions CShigh = Vehicle');
    ylabel('PExp Values');
    
    orthviews(fmri_data(stats{sigr}.weight_obj))
    
    % stats
    [tmain(sigr) pmain(sigr)] = ttest([nanmean(pexpvals_learning_A_CH{sigr}(:,1:2), 2); nanmean(pexpvals_learning_V_CH{sigr}(:,3:4), 2)], [nanmean(pexpvals_learning_A_CH{sigr}(:,3:4), 2); nanmean(pexpvals_learning_V_CH{sigr}(:,1:2), 2)]);
    
    input('Enter to continue')
    
end




    

    
%% use those 31 regions to test on generalization data? 
% load generalization data

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


%% apply pattern 

cmgen = [0 .9 .3; 0 .8 .2; 0 .7 .1; 0 .8 .6; 0 .7 .5; 0 .6 .4; 0 .7 .9; 0 .6 .8; 0 .5 .7; ...
         .9 0 .3; .8 0 .2; .7 0 .1; .8 0 .6; .7 0 .5; .6 0 .4; .7 0 .9; .6 0 .8; .5 0 .7];


for sigr = 1:numel(find(sign01))
    
    pexpvals_GCA1 = apply_mask(GCA1, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GCA2 = apply_mask(GCA2, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GCA3 = apply_mask(GCA3, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GPA1 = apply_mask(GPA1, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GPA2 = apply_mask(GPA2, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GPA3 = apply_mask(GPA3, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GWA1 = apply_mask(GWA1, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GWA2 = apply_mask(GWA2, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GWA3 = apply_mask(GWA3, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 

    pexpvals_GCV1 = apply_mask(GCV1, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GCV2 = apply_mask(GCV2, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GCV3 = apply_mask(GCV3, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GPV1 = apply_mask(GPV1, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GPV2 = apply_mask(GPV2, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GPV3 = apply_mask(GPV3, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GWV1 = apply_mask(GWV1, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GWV2 = apply_mask(GWV2, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 
    pexpvals_GWV3 = apply_mask(GWV3, stats{sigr}.weight_obj, 'pattern_expression', 'ignore_mising'); 

    pexpvals_gen{sigr} = [pexpvals_GCA1 pexpvals_GCA2 pexpvals_GCA3 ...
                    pexpvals_GPA1 pexpvals_GPA2 pexpvals_GPA3 ...
                    pexpvals_GWA1 pexpvals_GWA2 pexpvals_GWA3 ...
                    pexpvals_GCV1 pexpvals_GCV2 pexpvals_GCV3 ...
                    pexpvals_GPV1 pexpvals_GPV2 pexpvals_GPV3 ...
                    pexpvals_GWV1 pexpvals_GWV2 pexpvals_GWV3];


    clear pexpvals_G*
    
    create_figure('Categories_AVPatternexp GENERAL');
    barplot_colored(pexpvals_gen{sigr}, 'Pattern Expression', [], 'y',(-1:1)+.5,'nofig', 'noind', 'within');
    set(gca, 'XTick',[1:18], 'XTickLabel',{'', 'ANI cartoons', '', '', 'ANI photos', '', '', 'ANI words', '', '', 'VEH cartoons', '', '', 'VEH photos', '','', 'VEH words', ''});
    cm = [1 .5 0; .5 0 1; 0 .5 .2; 0 .1 .8];
    colormap(cmgen)
    xlabel('Session Conditions');
    ylabel('PExp Values');
    
    orthviews(fmri_data(stats{sigr}.weight_obj))
    
    input('Enter to continue')
    
    close all

end




%% do some stats on generalization classification performance

for sigr = find(t)

%     [t(sigr) p(sigr)] = ttest(nanmean(pexpvals_gen{sigr}(:,1:9), 2), nanmean(pexpvals_gen{sigr}(:,10:18), 2));
    
    create_figure('Categories_AVPatternexp GENERAL');
    barplot_colored(pexpvals_gen{sigr}, 'Pattern Expression', [], 'y',(-1:1)+.5,'nofig', 'noind', 'within');
    set(gca, 'XTick',[1:18], 'XTickLabel',{'', 'ANI cartoons', '', '', 'ANI photos', '', '', 'ANI words', '', '', 'VEH cartoons', '', '', 'VEH photos', '','', 'VEH words', ''});
    cm = [1 .5 0; .5 0 1; 0 .5 .2; 0 .1 .8];
    colormap(cmgen)
    xlabel('Session Conditions');
    ylabel('PExp Values');
    
    orthviews(fmri_data(stats{sigr}.weight_obj))
    
    input('Enter to continue')
    
end

