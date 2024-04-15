% Specify the paths
projectDir = '/home/dsutterlin/projects/genPain/';
canlabCorePath = fullfile(projectDir,'/canlab/CanlabCore');
mediationToolboxPath = fullfile(projectDir, '/canlab/MediationToolbox');
tutorialPath = fullfile(projectDir,'canlab/CANlab_help_examples/canlab_mixed_effects_matlab_demo1');

% Add paths to MATLAB
addpath(canlabCorePath);
addpath(mediationToolboxPath);
addpath(tutorialPath);

% Load the study_behav_dat_sid_14_164.csv file
data = readtable(fullfile(tutorialPath, 'study_behav_dat_sid_14_164.csv'));

% Now, you can run your tutorial code
%canlab_mixed_model_example;
behav_dat = readtable('study_behav_dat_sid_14_164.csv');

behav_dat(1:5, :)
behav_dat(behav_dat.sid == 113,:) = []
behav_dat = renamevars(behav_dat, 'prodicaine', 'placebo');

behav_dat(1:5,:)
heat_dat = behav_dat(behav_dat.heat == 1,:); % Use heat trials only

m = fitlme(heat_dat, 'Yint ~ stimLvl*placebo + (stimLvl*placebo | sid)',...
    'FitMethod','REML')
m.Coefficients

anova(m,'dfmethod','satterthwaite')


[Pa, Fa, DF1a, DF2a] = coefTest(m,[0,1,0,0],0,'dfmethod','satterthwaite');
disp(array2table([Pa, Fa, DF1a, DF2a], 'VariableNames', {'pValue', 'F', 'DF1', 'DF2'}, 'RowNames', {'stimLvl'}));


[Pa, Fa, DF1a, DF2a] = coefTest(m,[0,1,1,0],0,'dfmethod','satterthwaite');
disp('---stim and placebo---')
disp(array2table([Pa, Fa, DF1a, DF2a], 'VariableNames', {'pValue', 'F', 'DF1', 'DF2'}, 'RowNames', {'placebo = (-1)*stimLvl'}));

[Pa, Fa, DF1a, DF2a] = coefTest(m,[0,0,0,1],0,'dfmethod','satterthwaite');
disp('---Interaction---')
disp(array2table([Pa, Fa, DF1a, DF2a], 'VariableNames', {'pValue', 'F', 'DF1', 'DF2'}, 'RowNames', {'Interaction'}));


t = array2table([Pa,Fa,DF1a,DF2a; Pb, Fb, DF1b, DF2b],'VariableNames',...
    {'pValue','F','DF1','DF2'},'RowNames',{'stimLvl','placebo = (-1)*stimLvl'});

%%%%%%%%%%%%%%%%%%%%%%%%%555
help LinearMixedModel
m.Rsquared
m.VariableInfo                          % all variables in input dataset and which are in model
m.Variables(1:5, :)                     % all variables in input dataset
X = designMatrix(m); X(1:5, :)          % Fixed-effects design matrix
Xg = designMatrix(m, 'Random'); whos('Xg')          % Random-effects design matrix


% This option,  'DummyVarCoding', 'effects', effects-codes, which we
% prefer:
m = fitlme(heat_dat,'Yint ~ stimLvl*placebo + (stimLvl*placebo | sid)',...
    'FitMethod','REML', 'DummyVarCoding', 'effects');
% help anova says:
%     To obtain tests for the Type III hypotheses, set the 'DummyVarCoding'
%     to 'effects' when fitting the model.

% MULTI_LEVEL MODEL

u = unique(heat_dat.sid);
[Y, X] = deal(cell(1, length(u)));

Y_name = 'Yint';
X_var_names = {'stimLvl' 'placebo'};

Y_var = heat_dat.(Y_name);

X_var = zeros(size(Y_var, 1), length(X_var_names));

for j = 1:length(X_var_names)
    X_var(:, j) = heat_dat.(X_var_names{j});
end

for i = 1:length(u)

    wh = heat_dat.sid == u(i);  % indicator for this subject

    Y{i} = Y_var(wh);

    for j = 1:length(X_var_names)

        X{i}(:, j) = X_var(wh, j);

    end

end

% X2 is the 2nd-level (between-person) design matrix.
% we don't have any between-person predictors, so X2 is empty here:
X2 = [];
X2_names = {'2nd-level Intercept (average within-person effects)'};

stats = glmfit_multilevel(Y, X, X2, 'names', X_var_names, 'beta_names', X2_names);

% The 'weighted' option uses Empirical Bayes weighting for a single-step
% precision reweighting. This is a simple version of what is acccomplished
% by fitlme in Matlab or LMER in R, or igls.m (Lindquist/Wager).

stats = glmfit_multilevel(Y, X, X2, 'names', X_var_names, 'beta_names', X2_names, 'weighted');

% The 'boot' option, with 'bootsamples', n, runs n boostrap samples at the 2nd level:
% Run with 1000 bootstrap samples initially, and more are added based on P-value tolerance to avoid instability with too-few samples

stats = glmfit_multilevel(Y, X, X2, 'names', X_var_names, 'beta_names', X2_names, 'weighted', 'boot');

%stats = glmfit_multilevel(Y, X, X2, 'names', X_var_names, 'beta_names', X2_names, 'weighted', 'boot', 'plots');

effectcolors = {[.5 .5 .5] [.7 .4 .4] [.4 .4 .7]};

create_figure('2nd level stats');
barplot_columns(stats.first_level.beta', 'names', stats.inputOptions.names, 'colors', effectcolors, 'nofigure');
ylabel(Y_name);
title('Naive individual within-person scores');
%drawnow, snapnow

X_placebo = X;
for i = 1:length(X), X_placebo{i} = X_placebo{i}(:, 2); end
create_figure('1st level effect');
line_plot_multisubject(X_placebo, Y);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'Control' 'placebo'});
%
indiv_scores = stats.Y_star(:, 3);
% or: indiv_scores = stats.first_level.beta(3, :)';

wh = indiv_scores > 0; % individuals with positive effects

poscolors = {[1 .4 .4]}; % custom_colors([.3 .3 .7], [.3 .5 1], sum(wh))';    % colors for those with pos effects
negcolors = {[.3 .3 1]}; % custom_colors([.7 .2 .2], [.7 .4 .4], sum(~wh))';    % colors for everyone else

colors = {};
colors(wh) = poscolors;
colors(~wh) = negcolors;

create_figure('1st level effect');
line_plot_multisubject(X_placebo, Y, 'colors', colors);
set(gca, 'XTick', [-.5 .5], 'XLim', [-.6 .6], 'XTickLabel', {'Control' 'placebo'});

X2 = cellfun(@nanmean, Y)';

stats = glmfit_multilevel(Y, X, X2, 'names', X_var_names, 'beta_names', {'Average_pain'}, 'weighted', 'boot')
