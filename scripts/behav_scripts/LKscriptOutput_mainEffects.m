%glmfit_multilevel from canlab toolbox for main effects of learning and
%gen. task


clc
close all
clear all

projectDir = '/home/dsutterlin/projects/genPain/';
dataDir = fullfile(projectDir, 'DATA/Behavioral/')
saveDir = fullfile(projectDir, 'results/behavioral')
addpath(genpath(fullfile(projectDir, 'Toolboxes/')));

cd (saveDir);
load finalPreproc_simGen_data.mat SCEBLmri_gendata
%load simrate.mat
load SCEBLmri_learndata_FINAL_N36.mat SCEBLmri_data
%load SCEBLmri_data_gen36.mat
gdata_df = readtable('SCEBLmri_Gendata_TxT_N36.csv');
ldata_df = readtable('SCEBLmri_Learndata_TxT_N36.csv');

% loading and rename learn and gen data for regression
ldata = SCEBLmri_data;
gdata = SCEBLmri_gendata;

subjects = filenames(fullfile(dataDir,'SCEBL_MRI2*'));
effectcolors = {[.5 .5 .5] [.7 .4 .4] [.4 .4 .7]};
u = 1:36; %sub id
%[Y, X] = deal(cell(1, length(u)));

%% LEARN
%% EXP ~ CUES
% Extract Learners form ! Learners
%=========
Y_name = 'exp';
X_names = {'cues'};
X = ldata.cues
X2_names = {'2nd-level Intercept (average within-person effects)'};
    
cue_exp_eff = glmfit_multilevel(ldata.(Y_name), X, [], 'names', X_names, 'beta_names', X2_names, 'weighted')
learn_beta = cue_exp_eff.first_level.beta(2,:)'; 
learners = find(learn_beta > median(learn_beta));
nonlearners = find(learn_beta <= median(learn_beta));
learners_idx = learn_beta > median(learn_beta);
%Fig
create_figure('2nd level stats');
barplot_columns(cue_exp_eff.first_level.beta', 'names', cue_exp_eff.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Within-person effect of learning cue on pain expectation');
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(X, ldata.(Y_name));
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

% plot learner vs non-learners slopes
create_figure('Learners cue effect on pain');
line_plot_multisubject(X(learners), ldata.pain(learners));
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

create_figure('Non-Learners cue effect on pain');
line_plot_multisubject(X(nonlearners), ldata.pain(nonlearners));
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});


%% PAIN ~ CUES of 49 deg (medium) trials
%======================================
%Apply filter on trials to keep only trials with medium temp.
for sub = 1:numel(subjects)
    mask = ldata.temp{sub} == 49;
    ldata.cues49deg{sub} = ldata.cues{sub}(mask);
    ldata.pain49deg{sub} = ldata.pain{sub}(mask);
end
Y_name = 'pain';
X_var_names = {'Cues from medium pain trials'};
X2_names = {'2nd-level Intercept (average within-person effects)'};
X2 = learn_beta
cuePain = glmfit_multilevel(ldata.pain49deg, ldata.cues49deg, [], 'names', X_var_names, 'beta_names', X2_names, 'weighted')
create_figure('2nd level stats');
barplot_columns(cuePain.first_level.beta', 'names', cuePain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow

% Colors gradient as a function of learning
min_beta = min(learn_beta);
max_beta = max(learn_beta);
color_scale = colormap(parula); % Using parula colormap for a smooth transition from light yellow to orange
normalized_beta = (learn_beta - min_beta) / (max_beta - min_beta);
color_indices = ceil(normalized_beta * (size(color_scale, 1) - 1)) + 1;
rgb_colors = cell(size(color_indices));
for i = 1:length(color_indices)
    rgb_colors{i} = color_scale(color_indices(i), :);
end
min_beta = round(min(learn_beta));
max_beta = round(max(learn_beta));

% Multiline with bar plots for mean CS -----------
create_figure('1st level effect');
line_plot_multisubject(ldata.cues, ldata.pain,'colors', rgb_colors, 'exclude_low_range_Y');
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});
hold on;
conditions = [-1, 1]; % Unique conditions
mean_pain49 = zeros(1, length(conditions)); % Preallocate mean pain array
filtered_49 = (ldata_df.temp == 49);
% Calculate mean pain scores for each condition
for i = 1:length(conditions)
    condition_indices = (ldata_df.cues == conditions(i)) & filtered_49;
    mean_pain49(i) = nanmean(ldata_df.pain(condition_indices));
end
% Add bar plots of the mean pain scores
bar_handle = bar(conditions, mean_pain49);
bar_handle.FaceColor = [0, 0.4470, 0.7410]; % Light blue color
bar_handle.FaceAlpha = 0.5; % Transparency
ylim([0 80]);
xlabel('Conditioned cue conditions');
ylabel('Pain score');
title('Conditioned cues effect on pain');

hLegend = legend('Individual Slopes');
hold on;
p = patch(NaN, NaN, [0, 0.4470, 0.7410], 'FaceAlpha', 0.5);
set(get(get(p, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'on');
colormap(parula);
cbar = colorbar('Location', 'eastoutside');
cbar.Ticks = linspace(0, 1, 10);
rounded_tick_labels = round(linspace(min_beta, max_beta, 10)); % Round the tick labels
cbar.TickLabels = num2cell(rounded_tick_labels);
cbar.Label.String = 'Explicit learning';
hold off;

% ! MAin effect sig.


% PAIN ~ CUES * EXP !! NEED fitlme()                    
%==================()
Y_name = 'pain';
X_names = {'cues', 'exp'};
% Adjust struct based on VIs
X_cuesExp = cell(1, numel(subjects));
for sub = 1:numel(subjects)
    v1 = ldata.(X_names{1}){sub};
    v2 = ldata.(X_names{2}){sub};
    X_cuesExp{sub} = horzcat(v1, v2);cd 
end

X2_names = {'2nd-level Intercept (average within-person effects)'};

cue_exp_eff = glmfit_multilevel(ldata.(Y_name), X_cuesExp, [], 'names', X_names, 'beta_names', X2_names, 'weighted')

% Extract single col from X
X_cues = X_cuesExp;
for i = 1:length(X_cuesExp), X_cues{i} = X_cues{i}(:, 2); end
X_exp = X_cuesExp;
for i = 1:length(X_cuesExp), X_exp{i} = X_exp{i}(:,1); end
%fig
create_figure('2nd level stats');
barplot_columns(cue_exp_eff.first_level.beta', 'names', cue_exp_eff.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(X_cues, ldata.(Y_name));
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

create_figure('1st level effect');
line_plot_multisubject(X_exp, ldata.(Y_name));
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});


% EXP ~ CUES from 49deg (medium pain) trials
%===========================================
Y_name = 'exp';
X_names = {'cues'};
X2_names = {'2nd-level Intercept (average within-person effects)'};

cuePain = glmfit_multilevel(ldata.(Y_name), ldata.cues, [], 'names', X_var_names, 'beta_names', X2_names, 'weighted')
create_figure('2nd level stats');
barplot_columns(cuePain.first_level.beta', 'names', cuePain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(ldata.cues, ldata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});


% PAIN ~ TEMP
Y_name = 'pain';
X_var_names = {'temp'};
X2_names = {'2nd-level Intercept (average within-person effects)'};

cuePain = glmfit_multilevel(ldata.pain, ldata.temp, [], 'names', X_var_names, 'beta_names', X2_names, 'weighted', 'plots')
create_figure('2nd level stats');
barplot_columns(cuePain.first_level.beta', 'names', cuePain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(ldata.cues, ldata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});


%% GENERALIZATION TASK
%====================

plot_correlation_matrix(gdata_df, 'names', gdata_df.Properties.VariableNames)

%% PAIN ~ CUES (cscat)
%============
Y_name = 'pain';
X_names = {'cscat'}
X2_names = {'2nd-level Intercept (average within-person effects)'};

genPain = glmfit_multilevel(gdata.pain, gdata.cscat, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
gen_beta = genPain.first_level.beta(2,:)'; 

create_figure('2nd level stats');
barplot_columns(genPain.first_level.beta', 'names', genPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Within-person cue effect on pain ratings');
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(gdata.cscat, gdata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});

% GENERALIZATION EFFECTS ~ Exp
% Create a scatter plot
figure;
scatter(learn_beta, gen_beta, "filled");
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Generalization effect in relation to explicit expectations');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, gen_beta);
annotation('textbox', [0.62, 0.15, 0.1, 0.1], 'String', sprintf('Pearson r = : %.2f', r), 'FontSize', 12, 'EdgeColor', 'none');
annotation('textbox', [0.62, 0.08, 0.1, 0.1], 'String', sprintf('p = : %.4f', p), 'FontSize', 12, 'EdgeColor', 'none');

% Learners only
Y_name = 'pain';
X_names = {'cscat'}
X2_names = {'2nd-level Intercept (average within-person effects)'};

genPain = glmfit_multilevel(gdata.pain(learners_idx), gdata.cscat(learners_idx), [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
gen_beta = genPain.first_level.beta(2,:)'; 
create_figure('2nd level stats');
barplot_columns(genPain.first_level.beta', 'names', genPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);

% ---------------
% Gen cues model with second level expectations
X_names = {'cscat'};
X2 = [learn_beta];
genPain = glmfit_multilevel(gdata.pain, gdata.cscat, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plot')
create_figure('2nd level stats');
barplot_columns(genPain.Y_star, 'names', {'Intercept', 'cscat'}, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Individual within-person scores');
drawnow, snapnow

% Multiline with bar plots for mean CS -----------
create_figure('1st level effect');
line_plot_multisubject(gdata.cscat, gdata.pain,'colors', rgb_colors, 'exclude_low_range_Y');
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'Low pain' 'High pain'});
hold on;
conditions = [-1, 1]; % Unique conditions
mean_painGen = zeros(1, length(conditions)); % Preallocate mean pain array
%filtered_49 = (gdata_df.temp == 49);
% Calculate mean pain scores for each condition
for i = 1:length(conditions)
    condition_indices = (gdata_df.cscat == conditions(i));
    mean_painGen(i) = nanmean(gdata_df.pain(condition_indices));
end
% Add bar plots of the mean pain scores
bar_handle = bar(conditions, mean_painGen);
bar_handle.FaceColor = [0, 0.4470, 0.7410]; % Light blue color
bar_handle.FaceAlpha = 0.5; % Transparency
ylim([0 80]);
xlabel('Generalization conditions');
ylabel('Pain Score');
title('Generalization cues effect on pain');
hLegend = legend('Individual Slopes');
hold on;
p = patch(NaN, NaN, [0, 0.4470, 0.7410], 'FaceAlpha', 0.5);
set(get(get(p, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'on');

colormap(parula);
cbar = colorbar('Location', 'eastoutside');
cbar.Ticks = linspace(0, 1, 10);
rounded_tick_labels = round(linspace(min_beta, max_beta, 10)); % Round the tick labels
cbar.TickLabels = num2cell(rounded_tick_labels);
cbar.Label.String = 'Explicit learning';
hold off;


%% PAIN ~ Similarity
%==============
Y_name = 'pain';
X_names = {'empirical_sim'} %,'modality'};
simPain = glmfit_multilevel(gdata.pain, gdata.empiricalSimr, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
sim_betas = simPain.first_level.beta(2,:)'; 

create_figure('2nd level stats');
barplot_columns(simPain.first_level.beta', 'names', simPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Within-person similarity on pain ratings');
drawnow, snapnow

figure;
scatter(learn_beta, sim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning and similarity betas imilarity');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, sim_betas);
text(0.5, max(sim_betas) - 0.1, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, max(sim_betas) - 0.11, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% Mean model similarity 
Y_name = 'pain';
X_names = {'mean_modelSimr'} %,'modality'};
simPain = glmfit_multilevel(gdata.pain, gdata.mean_modelSimr, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
meanSim_betas = simPain.first_level.beta(2,:)'; 

figure;
scatter(learn_beta, meanSim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning and mean theoretical similarity betas');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, meanSim_betas);
text(0.5, max(sim_betas) - 1, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, max(sim_betas) - 1.2, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% boulder LSA similarity
Y_name = 'pain';
X_names = {'boulderLSASimr'} %,'modality'};
simPain = glmfit_multilevel(gdata.pain, gdata.boulderLSASimr, [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
meanSim_betas = simPain.first_level.beta(2,:)'; 

figure;
scatter(learn_beta, meanSim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('LSA (BoulderU) conceptual similarity and expectation');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, meanSim_betas);
text(0.5, -2, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, max(sim_betas) - 3, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% Learners only!
%==============
Y_name = 'pain';
X_names = {'empirical_sim'} %,'modality'};
simPain = glmfit_multilevel(gdata.pain(learners_idx), gdata.empiricalSimr(learners_idx), [], 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
sim_betas = simPain.first_level.beta(2,:)'; 

create_figure('2nd level stats');
barplot_columns(simPain.first_level.beta', 'names', simPain.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);


%% PAIN ~ CUES (cscat) + simR from learners only
%=======================================
S
drawnow, snapnow
create_figure('1st level effect');
line_plot_multisubject(gdata.cscat, gdata.pain);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'CSLow' 'CSHigh'});



% PAIN ~ Similarity and second level expectation : Not sig!
X_names = {'empiricalSimr'}; %,'modality'};
X2_names = {'Learning beta'};
X2 = learn_beta
genPain = glmfit_multilevel(gdata.pain, gdata.boulderSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot', 'plots')
gen_beta = genPain.first_level.beta(2,:)'; 

%% PAIN ~  CSCAT + simmilarity (trials)
%=====================
Y_name = 'pain';
X_names = {'cscat', 'mean_modelSimr'}; %theotrialSim
for sub = 1:numel(subjects);
    v1 = gdata.(X_names{1}){sub};
    v2 = gdata.(X_names{2}){sub};
    X_cat_trials{sub} = horzcat(v1, v2);
end
X2_names = {'2d level'};
X2 = []
genPain = glmfit_multilevel(gdata.pain, X_cat_trials, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')
gen_beta = genPain.first_level.beta(2,:)'; 
sim_betas = genPain.first_level.beta(3,:)';
% scatter plot cue betas
figure;
scatter(learn_beta, gen_beta);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning and generalization betas, controlling for similarity');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, gen_beta);
text(0.5, 3, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, 2.5, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% scatter plot sim betas
figure;
scatter(learn_beta, sim_betas);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Learning and similarity betas, controlling for CS effect');
hold on;
lsline;
hold off;
[r, p] = corr(learn_beta, sim_betas);
text(0.5, 3, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, 2.5, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

% Second level expectatation
X2 = [learn_beta]
genPain = glmfit_multilevel(gdata.pain, X_cat_trials, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')


plot_correlation_matrix(genPain.Y_star, 'names', genPain.inputOptions.names);
drawnow, snapnow
figure; plotmatrix(genPain.Y_star);
drawnow, snapnow
% GENERALIZATION EFFECTS ~ Exp
% Create a scatter plot
figure;
scatter(learn_beta, gen_beta);
xlabel('Explicit Learning (beta)');
ylabel('Generalization (beta)');
title('Scatter Plot');
% Add regression line
hold on;
lsline;
hold off;
% Calculate correlation coefficient and p-value
[r, p] = corr(learn_beta, gen_beta)

% Display correlation coefficient and p-value on the plot
text(0.5, 3, sprintf('Correlation: %.2f', r), 'FontSize', 12);
text(0.5, 2.5, sprintf('P-value: %.4f', p), 'FontSize', 12, 'Color', 'red');

%% Second level learning
% Second level expectation on Pain~sim
X_names = {'empiricalSimr'};
X2 = [learn_beta];
genPain = glmfit_multilevel(gdata.pain, gdata.empiricalSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')

X_names = {'mean_modelSimr'};
X2 = [learn_beta];
genPain = glmfit_multilevel(gdata.pain, gdata.mean_modelSimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')

X_names = {'boulderLSASimr'};
X2 = [learn_beta];
genPain = glmfit_multilevel(gdata.pain, gdata.boulderLSASimr, X2, 'names', X_names, 'beta_names', X2_names, 'weighted','boot')


%% Correl between var
colnames = {'blocktrial', 'condition','cscat', 'pain','empiricalSimr', 'mean_modelSimr', 'boulderLSASimr'}
plot_correlation_matrix(gdata_df(:,colnames), 'names', colnames);
drawnow, snapnow

%% Mediation models
M = gdata.empiricalSimr

[paths, stats] = mediation(M,gdata.pain,gdata.cscat, 'plots','boot', 'verbose', 'bootsamples', 10000);
mediation_path_diagram(stats)

[paths, stats] = mediation(gdata.cscat(learners),gdata.pain(learners),M(learners) , 'plots','boot', 'verbose', 'bootsamples', 10000);
mediation_path_diagram(stats)

[r,p] = corr(gdata_df.cscat, gdata_df.empiricalSimr)
[r,p] = corr(gdata_df.cscat, gdata_df.mean_modelSimr)
[r,p] = corr(gdata_df.cscat, gdata_df.boulderLSASimr)
%correl CS and sim model?
%mediation

/home/dylan.sutterlinguindon/genPain/behav_scripts/LKscripts/analyze_main_SCEBL_fMRI.m
/home/dylan.sutterlinguindon/genPain/results 
/crnldata/socialhealth/projects/2024_PainGen/Behavioral/DATA

