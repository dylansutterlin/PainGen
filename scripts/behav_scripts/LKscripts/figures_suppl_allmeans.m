%% Figures for SCEBL1 manuscript

clc
clear all
close all


%% 
timevec = -2:0.02:20;
cm = [.2 .2 .2; .4 .4 .4; 0 0 .6; 0 0 .9; .7 0 0; 1 0 0; .5 0 .2; 1 0 .4];
mconditions = {'CLSL47'; 'CLSH47'; 'CLSL4748'; 'CLSH48'; 'CHSL48'; 'CHSH48'; 'CHSL49'; 'CHSH49'};

set(gcf, 'DefaultLineLineWidth', 1.5)
set(0, 'DefaultFigureColor', [1 1 1])

%% Exp 2 (SCEBL_MRI)

cd ('/Users/leoniekoban/Dropbox/work/BOULDER_PROJECTS/6_SCEBL1_fMRI/Behavioral/Results/')
load SCEBLmri_data_FINAL_N36.mat
load meanratings_N36_FINAL.mat
load ('/Users/leoniekoban/Dropbox/work/BOULDER_PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/nps_mean.mat')
load ('/Users/leoniekoban/Dropbox/work/BOULDER_PROJECTS/6_SCEBL1_fMRI/Physio/Results/GSR_means.mat')

betas_mean = betas_mean([1:5 7:9 11:end],:);
nps_mean = nps_mean([1:5 7:9 11:end],:);

%% figure
% ----- behavior -----
figure('Color', [1 1 1])

subplot(1,3,1)
errorbar(48:49, mean(painmeans_all(:,[1 3])), repmat(barplot_get_within_ste(painmeans_all), 1, 2), '--', 'Color', cm(3,:), 'LineWidth', 1.5); hold on
errorbar(48:49, mean(painmeans_all(:,[2 4])), repmat(barplot_get_within_ste(painmeans_all), 1, 2), '-', 'Color', cm(4,:), 'LineWidth', 1.5); hold on
errorbar(49:50, mean(painmeans_all(:,[5 7])), repmat(barplot_get_within_ste(painmeans_all), 1, 2), '--', 'Color', cm(5,:), 'LineWidth', 1.5); hold on
errorbar(49:50, mean(painmeans_all(:,[6 8])), repmat(barplot_get_within_ste(painmeans_all), 1, 2), '-', 'Color', cm(6,:), 'LineWidth', 1.5); hold on
xlim([47.5 50.5])
ylim([10 60])
set(gca, 'XTick', 48:1:50, 'XTickLabel',{'48', '49', '50'}, 'FontName', 'Arial', 'FontSize', 16);  xlabel('Temperature'); ylabel('Pain rating')
% title('Exp 2 - Pain ratings', 'FontName', 'Arial', 'FontSize', 20);

% ---- physio -------
subplot(1,3,2)
errorbar(48:49, nanmean(betas_mean(:,[1 3])), repmat(barplot_get_within_ste(betas_mean), 1, 2), '--', 'Color', cm(3,:), 'LineWidth', 1.5); hold on
errorbar(48:49, nanmean(betas_mean(:,[2 4])), repmat(barplot_get_within_ste(betas_mean), 1, 2), '-', 'Color', cm(4,:), 'LineWidth', 1.5); hold on
errorbar(49:50, nanmean(betas_mean(:,[5 7])), repmat(barplot_get_within_ste(betas_mean), 1, 2), '--', 'Color', cm(5,:), 'LineWidth', 1.5); hold on
errorbar(49:50, nanmean(betas_mean(:,[6 8])), repmat(barplot_get_within_ste(betas_mean), 1, 2), '-', 'Color', cm(6,:), 'LineWidth', 1.5); hold on
xlim([47.5 50.5])
ylim([0 1])
set(gca, 'XTick', 48:1:50, 'XTickLabel',{'48', '49', '50'}, 'FontName', 'Arial', 'FontSize', 16); xlabel('Temperature'); ylabel('SCR beta')
% title('Exp 1 - SCR', 'FontName', 'Arial', 'FontSize', 20);

% % ----- NPS -----
subplot(1,3,3)
errorbar(48:49, mean(nps_mean(:,[1 3])), repmat(barplot_get_within_ste(nps_mean), 1, 2), '--', 'Color', cm(3,:), 'LineWidth', 1.5); hold on
errorbar(48:49, mean(nps_mean(:,[2 4])), repmat(barplot_get_within_ste(nps_mean), 1, 2), '-', 'Color', cm(4,:), 'LineWidth', 1.5); hold on
errorbar(49:50, mean(nps_mean(:,[5 7])), repmat(barplot_get_within_ste(nps_mean), 1, 2), '--', 'Color', cm(5,:), 'LineWidth', 1.5); hold on
errorbar(49:50, mean(nps_mean(:,[6 8])), repmat(barplot_get_within_ste(nps_mean), 1, 2), '-', 'Color', cm(6,:), 'LineWidth', 1.5); hold on
xlim([47.5 50.5])
ylim([20 100])
set(gca, 'XTick', 48:1:50, 'XTickLabel',{'48', '49', '50'}, 'FontName', 'Arial', 'FontSize', 16); xlabel('Temperature'); ylabel('NPS response')
% title('Exp 2 - NPS response', 'FontName', 'Arial', 'FontSize', 20);



