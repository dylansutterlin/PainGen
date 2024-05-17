
clc;
close all;
clear all;


projectDir = '/home/dsutterlin/projects/genPain/';
toolboxdir = fullfile(projectDir, 'Toolboxes');
dataDir = fullfile(projectDir, 'DATA/Imaging/test_genBetas/model3_generalization_ST');
behavDir = fullfile(projectDir,'results','behavioral');
saveDir = fullfile(projectDir, 'results/imaging/CVpred');

% cRNL cluster paths
projectDir = '/crnldata/socialhealth/projects/2024_PainGen/';
toolboxdir = fullfile(projectDir, 'Toolboxes');
dataDir = '/crnldata/socialhealth/data/2024_SCEBL_GEN/Imaging_Data/analyses/model3_generalization_ST';
behavDir = fullfile(projectDir,'results','behavioral');
saveDir = fullfile(projectDir, '/results/imaging/cvPred');

mkdir(saveDir);


% Add Toolboxes to paths
addpath(genpath(fullfile(toolboxdir, 'RobustToolbox')));
addpath(genpath(fullfile(toolboxdir, 'Neuroimaging_Pattern_Masks')));
addpath(genpath(fullfile(toolboxdir, 'CanlabCore')));
addpath(genpath(fullfile(toolboxdir, 'MediationToolbox')));
addpath(genpath(fullfile(toolboxdir, 'spm12')));

% Except SPM/external for compatibility with Canlab
% toolbox/mediation (?)
spmFolders = dir(fullfile(toolboxdir, 'spm12'));
for i = 1:numel(spmFolders)
    folderName = spmFolders(i).name;
    if strcmp(folderName, 'external')
        path2rm = strsplit(genpath(fullfile(toolboxdir, 'spm12', folderName)), ':')
        for i=1:numel(path2rm)
            rmpath(char(path2rm(i)));
        end
    else
        continue
    end
end


% --------------
% Load varialbes
% --------------
load(fullfile(behavDir, 'finalPreproc_simGen_data.mat'));
gendata = SCEBLmri_gendata;

load(fullfile(behavDir, 'SCEBLmri_learndata_FINAL_N36.mat'))
learndata = SCEBLmri_data;

load(fullfile(projectDir, '/results/imaging/STmodel3_cues_vifsQC/subjectsVifs.mat'));
vifStruct = subjects;

% Functionnal images
allFiles = filenames(fullfile(dataDir, 'SCEBL*'));
subjectsDir = [];
% remove .zip folder from files to extract func files
for i = 1:numel(allFiles)
    filename = allFiles(i);
    if ~endsWith(filename, '.zip')
        % Add the filename to the subjects cell array
        subjectsDir = [subjectsDir; filename];
    end
end

% Extract subjects names
subjects = cell(1, numel(subjectsDir));
for p = 1:numel(subjectsDir)
    [~, subname] = fileparts(subjectsDir(p));
    subjects{p} = subname;
end 

% Initialize an empty cell array one for each subject
func_imgs = cell(1, numel(subjects)); % Initialize func_imgs as a cell array

% get func files in a cell/sub
for sub = 1:numel(subjects)
    contrastFiles = {}; % Initialize contrastFiles as an empty cell array
    subFiles = fullfile(dataDir, char(subjects{sub}));
    
    SPM = load(fullfile(subFiles, 'SPM.mat')); % Load SPM.mat for the subject
    table = struct2table(SPM.SPM.Vbeta);

    % Loop through each contrast
    for cont = 1:numel(table.descrip)
        % Check if the description of the contrast contains 'allPains_trials0'
        if contains(table.descrip{cont}, 'AllPains_trial0')
            % add the image path to contrastFiles
            add_img = table.fname{cont};
            img_path = fullfile(dataDir, subjects{sub}, add_img);
            contrastFiles = [contrastFiles, img_path]; % Append image path to contrastFiles
        end
    end
    %func_imgs{sub} = contrastFiles;
    func_imgs{sub} = char(contrastFiles'); % Assign contrastFiles to the corresponding subject in func_imgs
end

%% Prep varialbes, scale and compute 2d level moderator
% X
cscat = gendata.cscat(1:26);
empSim = gendata.empiricalSimr(1:26);
meanSim = gendata.mean_modelSimr(1:26);
lsaSim = gendata.boulderLSASimr(1:26);
pain = gendata.pain(1:26); % Y
cov = gendata.blocktrial(1:26);
data = struct();
data.subjects = subjects;
data.func = func_imgs; % M
data.pain = pain;

% Mean center vairables
for i=1:numel(func_imgs) % mean center each cells
    x = cscat{i};
    data.cscat{i} = (x-mean(x))/std(x);

    x = empSim{i};
    data.empSim{i} = (x-mean(x))/std(x);

    x = meanSim{i};
    data.meanSim{i} = (x-mean(x))/std(x);

    x = lsaSim{i};
    data.lsaSim{i} = (x-mean(x))/std(x);

    x = cov{i};
    data.cov{i} = (x-mean(x))/std(x);

    %y=pain{i};
    %data.pain{i} = (y-mean(y))/std(y);
end 

% L2M : Expectation (during learning) as moderator (expectation ~ learning cues)
cue_exp_eff = glmfit_multilevel(learndata.exp, learndata.cues, [], 'weighted');
learn_beta = cue_exp_eff.first_level.beta(2,:)';
l2m = cell(1, numel(learn_beta));
for m = 1:numel(learn_beta)
    l2m{m} = learn_beta(m);
end
% mean center
data.learnBetas = (learn_beta-mean(learn_beta))/std(learn_beta);
data.txtSubjects = gendata.sub_TxTgrp;

%% Remove outlier from X M and Y based on VIF   

for sub = 1:length(data.subjects)
    mask = vifStruct.outvifs{sub}; % Vector of 0 (keep) or 1 (rm)
    keep_mask = mask == 0; % keep ones and leaves the 0; (inverse of mask)
    
    data.txtSubjects{sub} = data.txtSubjects{sub}(keep_mask);
    data.cscat{sub} = data.cscat{sub}(keep_mask);
    data.empSim{sub} = data.empSim{sub}(keep_mask);
    data.meanSim{sub} = data.meanSim{sub}(keep_mask);
    data.lsaSim{sub} = data.lsaSim{sub}(keep_mask);
    data.pain{sub} = data.pain{sub}(keep_mask);
    data.cov{sub} = data.cov{sub}(keep_mask);
    data.func{sub} = data.func{sub}(keep_mask, :);
end

%% CV predictions 
% ----------------------------------------

mkdir(fullfile(saveDir, 'svm_cscat'));
cd(fullfile(saveDir, 'svm_cscat'));

dat = fmri_data(data.func);
dat.Y = cat(1, data.cscat{:}); % transform cell to vector
subject_id = cat(1, gendata.sub_TxTgrp{:});           

holdout_set = zeros(size(subject_id));      % holdout set membership for each image
n_subjects = length(unique(subject_id));
C = cvpartition(n_subjects, 'KFold', 5);
for i = 1:5
    teidx = test(C, i);                     % which subjects to leave out
    imgidx = ismember(subject_id, find(teidx)); % all images for these subjects
    holdout_set(imgidx) = i;
end

algoname = 'cv_svm'; % cross-validated penalized regression. Predict pain ratings

[cverr, stats, optout] = predict(dat,  'cv_svm', 'nfolds', 5, 'error_type', 'mse' );

save(pwd, 'cverr', 'stats', 'optout');

%% lasso-PCR with similarity as Y
mkdir(fullfile(saveDir, 'lassoPCR_empSim'));
cd(fullfile(saveDir, 'svm_cscat'));

dat.Y = cat(1, data.empSim{:});   

holdout_set = zeros(size(subject_id));      % holdout set membership for each image
n_subjects = length(unique(subject_id));
C = cvpartition(n_subjects, 'KFold', 5);
for i = 1:5
    teidx = test(C, i);                     % which subjects to leave out
    imgidx = ismember(subject_id, find(teidx)); % all images for these subjects
    holdout_set(imgidx) = i;
end

algoname = 'cv_lassopcr'; % cross-validated penalized regression. Predict pain ratings
[cverr, stats, optout] = predict(dat,  algoname, 'nfolds', 5, 'bootweights', 'bootsamples', 100);


%% Figures
create_figure('ROC')
ROC = roc_plot(stats.yfit, logical(stats.Y > 0), 'color', 'r');

%[m,aprime,corr,far,missclass, mprop, statstable] = confusion_matrix(stats.Y, stats.yfit);
%mytable = array2table(mprop, 'VariableNames', {'Pred Cog' 'Pred Pain'}, 'RowNames', {'Actual Cog' 'Actual Pain'});
%disp(mytable)

dat_to_plot = {stats.dist_from_hyperplane_xval(stats.Y > 0) stats.dist_from_hyperplane_xval(stats.Y < 0)};
create_figure('SVM distance scores by class');
barplot_columns(dat_to_plot, 'title', 'Training data SVM scores', 'nofigure', 'colors', {[.7 .2 .2] [.5 .3 .7]});

create_figure('weight map'); axis off;
montage(stats.weight_obj);
