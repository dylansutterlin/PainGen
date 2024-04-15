
% addpath('/Users/Leonie/Documents/Tools/scripts/General scripts')

clc
close all
clear all

projectDir = '/home/dsutterlin/projects/genPain';
dataDir = fullfile(projectDir, 'DATA/Behavioral/')
saveDir = fullfile(projectDir, 'results/behavioral')

%load ('C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI/Behavioral/Results/SCEBLmri_data_FINAL_N36.mat');
load (fullfile(saveDir, 'SCEBLmri_learndata_FINAL_N36.mat'))

%load ('C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI/Behavioral/Results/simmat.mat');
load (fullfile(saveDir, 'simmat.mat')) 

cd (dataDir);

subjects = {'SCEBL_MRI201';'SCEBL_MRI202';'SCEBL_MRI203';'SCEBL_MRI204';'SCEBL_MRI205';...
                           'SCEBL_MRI207';'SCEBL_MRI208';'SCEBL_MRI209';               ...
            'SCEBL_MRI211';'SCEBL_MRI212';'SCEBL_MRI213';'SCEBL_MRI214';'SCEBL_MRI215';...
            'SCEBL_MRI216';'SCEBL_MRI217';'SCEBL_MRI218';...
            'SCEBL_MRI219';'SCEBL_MRI220';'SCEBL_MRI221';'SCEBL_MRI222';'SCEBL_MRI223';...
            'SCEBL_MRI224';'SCEBL_MRI225';'SCEBL_MRI226';'SCEBL_MRI227';'SCEBL_MRI228';...
            'SCEBL_MRI229';'SCEBL_MRI230';'SCEBL_MRI231';'SCEBL_MRI232';'SCEBL_MRI233';...
            'SCEBL_MRI234';'SCEBL_MRI235';'SCEBL_MRI236';'SCEBL_MRI237';'SCEBL_MRI238'};
        
% bad/missing data from 206 and 210; 


% whether A(1) or V(2) is CShigh
AV_h = [2, 1, 2, 1, 2,...
           2, 1, 2,   ...
        1, 2, 1, 2, 1, ...
        2, 1, 2, ...
        2, 1, 2, 1, 2,...
        1, 2, 1, 2, 1,...
        1, 2, 1, 2, 1,...
        2, 1, 2, 1, 1];

    
Ahigh = find(AV_h == 1);
Vhigh = find(AV_h == 2);

cue_exp_eff = glmfit_multilevel(SCEBLmri_data.exp, SCEBLmri_data.cues, []);
learn_beta = cue_exp_eff.first_level.beta(2,:)';
        
learners = find(learn_beta > median(learn_beta));
nonlearners = find(learn_beta <= median(learn_beta));

gencons = {'GCA1'; 'GCA2'; 'GCA3'; 'GCV1'; 'GCV2'; 'GCV3';...
           'GPA1'; 'GPA2'; 'GPA3'; 'GPV1'; 'GPV2'; 'GPV3'; ...
           'GWA1'; 'GWA2'; 'GWA3'; 'GWV1'; 'GWV2'; 'GWV3'};

genconslegend = {'GCA1'; 'GCA2'; 'GCA3'; 'GPA1'; 'GPA2'; 'GPA3'; 'GWA1'; 'GWA2'; 'GWA3'; ...
                 'GCV1'; 'GCV2'; 'GCV3'; 'GPV1'; 'GPV2'; 'GPV3'; 'GWV1'; 'GWV2'; 'GWV3'};

       
       
% depending on how many runs subjects did, indexes in excel file (normal is 4 runs)
ind_4runs = [1:18 20:37 39:56 58:75];
 
% colormap cmgen
cmgen = [0 .9 .3; 0 .8 .2; 0 .7 .1; 0 .8 .6; 0 .7 .5; 0 .6 .4; 0 .7 .9; 0 .6 .8; 0 .5 .7; ...
             .9 0 .3; .8 0 .2; .7 0 .1; .8 0 .6; .7 0 .5; .6 0 .4; .7 0 .9; .6 0 .8; .5 0 .7];     

for sub = 1:numel(subjects)

   % put current subject folder here:
   subject = (subjects{sub});

   %cd(['C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI/Behavioral/DATA/', subject, '/3_Main'])
   cd (fullfile(dataDir, subject, '/3_Main'))

   lf = filenames('L*.xlsx', 'char');
   SCEBLmri_gendata.sub_TxTgrp{sub} = sub*ones(72, 1); %added var to account for trial subj's id
   SCEBLmri_gendata.subjname{sub} = (subjects{sub});
   SCEBLmri_gendata.cb{sub} = xlsread(lf, [lf(1:4), '.txt'], 'D3');
   [dummy, SCEBLmri_gendata.lowcue{sub}] = xlsread(lf, [lf(1:4), '.txt'], 'N3');
   [dummy, SCEBLmri_gendata.highcue{sub}] = xlsread(lf, [lf(1:4), '.txt'], 'M3');
   
   %cd(['C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI/Behavioral/DATA/', subject, '/4_Gen'])
   cd (fullfile(dataDir, subject, '/4_Gen'))

   f = filenames('G*.xlsx', 'char');
   
   SCEBLmri_gendata.occid{sub} = xlsread(f, [f(1:4), '.txt'], 'B3');
   
   ind_4runs = [1:18 20:37 39:56 58:75];

   dummy = xlsread(f, [f(1:4), '.txt'], 'AX3:AX78'); 
   SCEBLmri_gendata.blocktrial{sub} = dummy(ind_4runs);

   dummy = xlsread(f, [f(1:4), '.txt'], 'AV3:AV78'); 
   SCEBLmri_gendata.condition{sub} = dummy(ind_4runs);
   
   SCEBLmri_gendata.AVcat{sub}(SCEBLmri_gendata.condition{sub} == 1 |...
                               SCEBLmri_gendata.condition{sub} == 2 |... 
                               SCEBLmri_gendata.condition{sub} == 3 |...
                               SCEBLmri_gendata.condition{sub} == 7 |...
                               SCEBLmri_gendata.condition{sub} == 8 |...
                               SCEBLmri_gendata.condition{sub} == 9 |...
                               SCEBLmri_gendata.condition{sub} == 13 |...
                               SCEBLmri_gendata.condition{sub} == 14 |...
                               SCEBLmri_gendata.condition{sub} == 15 , 1) = -1;
                           
   SCEBLmri_gendata.AVcat{sub}(SCEBLmri_gendata.condition{sub} == 4 |...
                               SCEBLmri_gendata.condition{sub} == 5 |... 
                               SCEBLmri_gendata.condition{sub} == 6 |...
                               SCEBLmri_gendata.condition{sub} == 10 |...
                               SCEBLmri_gendata.condition{sub} == 11 |...
                               SCEBLmri_gendata.condition{sub} == 12 |...
                               SCEBLmri_gendata.condition{sub} == 16 |...
                               SCEBLmri_gendata.condition{sub} == 17 |...
                               SCEBLmri_gendata.condition{sub} == 18 , 1) = 1;
   
   SCEBLmri_gendata.modality{sub}(SCEBLmri_gendata.condition{sub} < 7, 1) = 1;
   SCEBLmri_gendata.modality{sub}((SCEBLmri_gendata.condition{sub} > 6 & SCEBLmri_gendata.condition{sub} < 13), 1) = 2;
   SCEBLmri_gendata.modality{sub}(SCEBLmri_gendata.condition{sub} > 12, 1) = 3;
   
   if AV_h(sub) == 1 % if Animal is CShigh inverse AVcat
       SCEBLmri_gendata.cscat{sub}(SCEBLmri_gendata.AVcat{sub}== -1, 1) = 1;
       SCEBLmri_gendata.cscat{sub}(SCEBLmri_gendata.AVcat{sub}== 1, 1) = -1;
       
   elseif AV_h(sub) == 2
       SCEBLmri_gendata.cscat{sub}(SCEBLmri_gendata.AVcat{sub}== -1, 1) = -1;
       SCEBLmri_gendata.cscat{sub}(SCEBLmri_gendata.AVcat{sub}== 1, 1) = 1;
       
   end
   
   [dummy2, dummy] = xlsread(f, [f(1:4), '.txt'], 'U3:U78'); 
   SCEBLmri_gendata.conditionname{sub} = dummy(ind_4runs);

   SCEBLmri_gendata.sim_lc{sub} = simmat{sub}(1, SCEBLmri_gendata.condition{sub})'; % get sim rating for the condition of each trial
   SCEBLmri_gendata.sim_hc{sub} = simmat{sub}(2, SCEBLmri_gendata.condition{sub})';
   SCEBLmri_gendata.sim_diff{sub} = SCEBLmri_gendata.sim_hc{sub} - SCEBLmri_gendata.sim_lc{sub};

   dummy = xlsread(f, [f(1:4), '.txt'], 'AI3:AI78'); 
   SCEBLmri_gendata.RTpain{sub} = dummy(ind_4runs);

   dummy = xlsread(f, [f(1:4), '.txt'], 'AD3:AD78'); 
   SCEBLmri_gendata.pain{sub} = dummy(ind_4runs);

   
   %% figure
   figure('Color', [1 1 1], 'Name', subject)
   
   title(subject);

   subplot(3,2,1); 
   plot(1:18, SCEBLmri_gendata.pain{sub}(1:18)); hold on
   plot(19:36, SCEBLmri_gendata.pain{sub}(19:36)); hold on
   plot(37:54, SCEBLmri_gendata.pain{sub}(37:54)); hold on
   plot(55:72, SCEBLmri_gendata.pain{sub}(55:72)); hold on
   xlim([0 73])
     
   for con = 1:numel(gencons)
       pain.(gencons{con})(sub,:) = SCEBLmri_gendata.pain{sub}(strcmp(SCEBLmri_gendata.conditionname{sub}, gencons{con}))';
   end
%% 
   subplot(3,2,2)
   title(['CSlow = ', SCEBLmri_gendata.lowcue{sub}, '   CShigh = ', SCEBLmri_gendata.highcue{sub}])

   barplot_colored([pain.GCA1(sub,:); pain.GCA2(sub,:); pain.GCA3(sub,:);...
                    pain.GPA1(sub,:); pain.GPA2(sub,:); pain.GPA3(sub,:);...
                    pain.GWA1(sub,:); pain.GWA2(sub,:); pain.GWA3(sub,:);...
                    pain.GCV1(sub,:); pain.GCV2(sub,:); pain.GCV3(sub,:);...
                    pain.GPV1(sub,:); pain.GPV2(sub,:); pain.GPV3(sub,:);...
                    pain.GWV1(sub,:); pain.GWV2(sub,:); pain.GWV3(sub,:)]'); hold on

   set(gca, 'XTick', [1:18], 'XTickLabel', genconslegend) 

   plot(simmat{sub}(1,:), 'k--', 'LineWidth', 2);hold on;
   plot(simmat{sub}(2,:), 'm', 'LineWidth', 2); hold on;

   
   title([subject, ' ---- ', 'CSlow = ', SCEBLmri_gendata.lowcue{sub}, '   CShigh = ', SCEBLmri_gendata.highcue{sub}]);
   
   %%
   
   animalmean(sub,:) = nanmean([pain.GCA1(sub,:); pain.GCA2(sub,:); pain.GCA3(sub,:);...
                             pain.GPA1(sub,:); pain.GPA2(sub,:); pain.GPA3(sub,:);...
                             pain.GWA1(sub,:); pain.GWA2(sub,:); pain.GWA3(sub,:)]');
   vehiclemean(sub,:) = nanmean([pain.GCV1(sub,:); pain.GCV2(sub,:); pain.GCV3(sub,:);...
                              pain.GPV1(sub,:); pain.GPV2(sub,:); pain.GPV3(sub,:);...
                              pain.GWV1(sub,:); pain.GWV2(sub,:); pain.GWV3(sub,:)]');
   
   if AV_h(sub) == 1
       CShighCat_mean(sub,:) = animalmean(sub,:);
       CSlowCat_mean(sub,:) = vehiclemean(sub,:);
   elseif AV_h(sub) == 2
       CShighCat_mean(sub,:) = vehiclemean(sub,:);
       CSlowCat_mean(sub,:) = animalmean(sub,:);    
   end
   
   subplot(3,2,3)
   barplot_colored([animalmean(sub,:);  vehiclemean(sub,:)]');

   
   genpainmeans{sub} = nanmean([pain.GCA1(sub,:); pain.GCA2(sub,:); pain.GCA3(sub,:);...
                                        pain.GPA1(sub,:); pain.GPA2(sub,:); pain.GPA3(sub,:);...
                                        pain.GWA1(sub,:); pain.GWA2(sub,:); pain.GWA3(sub,:);...
                                        pain.GCV1(sub,:); pain.GCV2(sub,:); pain.GCV3(sub,:);...
                                        pain.GPV1(sub,:); pain.GPV2(sub,:); pain.GPV3(sub,:);...
                                        pain.GWV1(sub,:); pain.GWV2(sub,:); pain.GWV3(sub,:)],2);
   sim_lc{sub} = simmat{sub}(1,:)';
   sim_hc{sub} = simmat{sub}(2,:)';                                     

   subplot(3,2,4)
   imagesc(simmat{sub});
   colorbar
   colormap(cmgen)

   subplot(3,2,5)
   scatter(sim_lc{sub}, genpainmeans{sub})
   xlabel('Sim to CSlow'); ylabel('Pain rating')

   subplot(3,2,6)
   scatter(sim_hc{sub}, genpainmeans{sub})
   xlabel('Sim to CShigh'); ylabel('Pain rating')
    
   close all
   
end

cd(saveDir)
save SCEBLmri_data_gen36.mat SCEBLmri_gendata

% =======================================================================
% =======================================================================
% =======================================================================
%% 

figure('Color', [1 1 1])
subplot(5,1,1)
barplot_colored([animalmean  vehiclemean], 'within')
title('all subjects', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 16)
set(gca, 'XTick', [1:18], 'XTickLabel', gencons, 'FontName', 'Arial', 'FontSize', 12) 

subplot(5,1,2)
barplot_colored([CSlowCat_mean  CShighCat_mean], 'within')
title('all subjects', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 16)
set(gca, 'XTick', [5 14], 'XTickLabel', {'CS High Cat' 'CS Low Cat'}, 'FontName', 'Arial', 'FontSize', 12) 

subplot(5,2,5)
barplot_colored([mean(animalmean(Vhigh,:), 2)  mean(vehiclemean(Vhigh,:),2)], 'within'); hold on
line([1 2], [mean(animalmean(Vhigh,:), 2)  mean(vehiclemean(Vhigh,:),2)])
title('Animal CSlow, Vehicle CShigh', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 14)
set(gca, 'XTick', [1:2], 'XTickLabel', {'Animals' 'Vehicles'}, 'FontName', 'Arial', 'FontSize', 14) 

subplot(5,2,6)
imagesc(nanmean(simmat2(:,:,Vhigh), 3))
title('Animal CSlow, Vehicle CShigh, Similarity matrix', 'FontName', 'Arial', 'FontSize', 16); colorbar
set(gca, 'YTick', [1:2], 'YTickLabel', {'CSlow' 'CShigh'}, 'FontName', 'Arial', 'FontSize', 14) 

subplot(5,2,7)
barplot_colored([mean(animalmean(Ahigh,:), 2)  mean(vehiclemean(Ahigh,:),2)], 'within'); hold on
line([1 2], [mean(animalmean(Ahigh,:), 2)  mean(vehiclemean(Ahigh,:),2)])
title('Animal CShigh, Vehicle CSlow', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 14)
set(gca, 'XTick', [1:2], 'XTickLabel', {'Animals' 'Vehicles'}, 'FontName', 'Arial', 'FontSize', 14) 

subplot(5,2,8)
imagesc(nanmean(simmat2(:,:,Ahigh), 3))
title('Animal CShigh, Vehicle CSlow, Similarity matrix', 'FontName', 'Arial', 'FontSize', 16); colorbar
set(gca, 'YTick', [1:2], 'YTickLabel', {'CSlow' 'CShigh'}, 'FontName', 'Arial', 'FontSize', 14) 

subplot(5,3,13)
barplot_colored([mean(CSlowCat_mean, 2)  mean(CShighCat_mean, 2)], 'within')
title('By CS category', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 14)
set(gca, 'XTick', [1:2], 'XTickLabel', {'CSlow Category' 'CShigh Category'}, 'FontName', 'Arial', 'FontSize', 14) 

subplot(5,3,14)
dummy = [mean(CSlowCat_mean, 2)  mean(CShighCat_mean, 2)];
barplot_colored(dummy(nonlearners,:), 'within')
title('Non-Learners by CS category', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 14)
set(gca, 'XTick', [1:2], 'XTickLabel', {'CSlow Category' 'CShigh Category'}, 'FontName', 'Arial', 'FontSize', 14) 

subplot(5,3,15)
barplot_colored(dummy(learners,:), 'within')
title('Learners by CS category', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 14)
set(gca, 'XTick', [1:2], 'XTickLabel', {'CSlow Category' 'CShigh Category'}, 'FontName', 'Arial', 'FontSize', 14) 

% [t p] = ttest([mean(animalmean(Ahigh,:), 2); mean(vehiclemean(Vhigh,:), 2)], [mean(vehiclemean(Ahigh,:),2);  mean(animalmean(Vhigh,:),2)]);
colormap(cmgen); 

figure
barplot_colored([mean(animalmean, 2)  mean(vehiclemean, 2)], 'within'); hold on
line([1 2], [mean(animalmean, 2)  mean(vehiclemean,2)])
title('All subjects', 'FontName', 'Arial', 'FontSize', 16); ylabel('Pain rating', 'FontName', 'Arial', 'FontSize', 14)
set(gca, 'XTick', [1:2], 'XTickLabel', {'Animals' 'Vehicles'}, 'FontName', 'Arial', 'FontSize', 14) 

[h p] = ttest(mean(animalmean, 2), mean(vehiclemean, 2)); % no sign. diff between animals and vehicles -- ooufff :)


% %% 
% figure('Color', [1 1 1])
% subplot(1,2,1)
% line_plot_multisubject(sim_lc, genpainmeans)
% xlabel('Sim to CSlow'); ylabel('Pain rating')
% 
% subplot(1,2,2)
% line_plot_multisubject(sim_hc, genpainmeans)
% xlabel('Sim to CShigh'); ylabel('Pain rating')
% 
% glmfit_multilevel(sim_hc, genpainmeans, learn_beta);
% 
% glmfit_multilevel(sim_hc, genpainmeans, AV_h'); % not significant, but also makes overall effect n.s.
% 
% % [r p] = corr(AV_h', learn_beta)  % n.s. luckily
% 
% 
% figure('Color', [1 1 1])
% subplot(2,2,1)
% line_plot_multisubject(sim_lc(nonlearners), genpainmeans(nonlearners))
% xlabel('NON-LEARNERS Sim to CSlow'); ylabel('Pain rating')
% 
% subplot(2,2,2)
% line_plot_multisubject(sim_hc(nonlearners), genpainmeans(nonlearners))
% xlabel('NON-LEARNERS Sim to CShigh'); ylabel('Pain rating')
% 
% subplot(2,2,3)
% line_plot_multisubject(sim_lc(learners), genpainmeans(learners))
% xlabel('LEARNERS Sim to CSlow'); ylabel('Pain rating')
% 
% subplot(2,2,4)
% line_plot_multisubject(sim_hc(learners), genpainmeans(learners))
% xlabel('LEARNERS Sim to CShigh'); ylabel('Pain rating')
% 
% 
% % %% load NPS
% 
% %load ('C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI/Behavioral/Results/pain_exp_gen36.mat');
% 
% SCEBLmri_gendata.nps = pain_exp;
% 
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.nps, []);
% figure; line_plot_multisubject(SCEBLmri_gendata.nps, SCEBLmri_gendata.pain);
% 
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.cscat, []);
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.modality, []);
% 
% 
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.sim_lc, []);
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.sim_hc, []);
% glmfit_multilevel(SCEBLmri_gendata.pain(learners), SCEBLmri_gendata.sim_hc(learners), []); % only learners
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.sim_diff, []);
% 
% glmfit_multilevel(SCEBLmri_gendata.nps, SCEBLmri_gendata.sim_lc, []);
% glmfit_multilevel(SCEBLmri_gendata.nps, SCEBLmri_gendata.sim_hc, []);
% glmfit_multilevel(SCEBLmri_gendata.nps, SCEBLmri_gendata.sim_diff, []);
% 
% % with 2nd level predictor
% 
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.cscat, learn_beta);
% 
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.nps, learn_beta);
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.sim_lc, learn_beta);
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.sim_hc, learn_beta);
% glmfit_multilevel(SCEBLmri_gendata.pain, SCEBLmri_gendata.sim_diff, learn_beta);
% 
% glmfit_multilevel(SCEBLmri_gendata.nps, SCEBLmri_gendata.sim_lc, learn_beta);
% glmfit_multilevel(SCEBLmri_gendata.nps, SCEBLmri_gendata.sim_hc, learn_beta);
% glmfit_multilevel(SCEBLmri_gendata.nps, SCEBLmri_gendata.sim_diff, learn_beta);
% 
% figure; barplot_colored(cell2mat(pain_exp)');
% 
% %cd C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI/Behavioral/Results
% 
% save SCEBLmri_data_gen36.mat SCEBLmri_gendata
% 
% %% remove high vif trials from single trial NPS data
% 
% % for s = 1:36
% %     SCEBLmri_gendata.high_vif_trials_idx{1,s} = find(painvifs{s} > 2.5);
% %     pain_exp{s}(SCEBLmri_gendata.high_vif_trials_idx{1,s}) = NaN;
% % end
% 
% nps_mat = cell2mat(pain_exp);
% 
% for subj = 1:36
% 
%     nps_mean(subj,1) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 1), subj));
%     nps_mean(subj,2) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 2), subj));
%     nps_mean(subj,3) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 3), subj));
%     nps_mean(subj,4) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 7), subj));
%     nps_mean(subj,5) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 8), subj));
%     nps_mean(subj,6) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 9), subj));
%     nps_mean(subj,7) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 13), subj));
%     nps_mean(subj,8) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 14), subj));
%     nps_mean(subj,9) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 15), subj));
%     nps_mean(subj,10) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 4), subj));
%     nps_mean(subj,11) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 5), subj));
%     nps_mean(subj,12) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 6), subj));
%     nps_mean(subj,13) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 10), subj));
%     nps_mean(subj,14) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 11), subj));
%     nps_mean(subj,15) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 12), subj));
%     nps_mean(subj,16) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 16), subj));
%     nps_mean(subj,17) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 17), subj));
%     nps_mean(subj,18) = nanmean(nps_mat((SCEBLmri_gendata.condition{subj} == 18), subj));
% 
% end
% 
% figure('Color', [1 1 1]); barplot_colored(nps_mean, 'within'); hold all;
% for cc = 1:18
%     plot(cc, nps_mean(:,cc), 'ko'); hold on;
% end
% 
% set(gca, 'XTick', [1:18], 'XTickLabel', genconslegend);
% 
% figure('Color', [1 1 1], 'Name', 'NPS responses')
% for subj = 1:36
%     subplot(6,6, subj)
%     bar(nps_mean(subj,:))
%     title(subjects{subj})
% end
% 
% %% by CS category
% 
% for sub = 1:36
%    if AV_h(sub) == 1
%        CShighCat_mean_NPS(sub,:) = nps_mean(sub,1:9);
%        CSlowCat_mean_NPS(sub,:) = nps_mean(sub,10:18);
%    elseif AV_h(sub) == 2
%        CShighCat_mean_NPS(sub,:) = nps_mean(sub,10:18);
%        CSlowCat_mean_NPS(sub,:) = nps_mean(sub,1:9);    
%    end
% end
% 
% figure('Color', [1 1 1]); barplot_colored([CSlowCat_mean_NPS CShighCat_mean_NPS], 'within'); hold all;
% bycatNPS = [CSlowCat_mean_NPS CShighCat_mean_NPS];
% for cc = 1:18
%     plot(cc, bycatNPS(:,cc), 'ko'); hold on;
% end
% 
% figure('Color', [1 1 1], 'Name', 'NPS responses')
% for subj = 1:36
%     subplot(6,6, subj)
%     bar([CSlowCat_mean_NPS(subj,:) CShighCat_mean_NPS(subj,:)])
%     title(subjects{subj})
% end



