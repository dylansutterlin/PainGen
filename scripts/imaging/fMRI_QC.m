% QC check

clc;
close all;
clear all;

%% Path prep
% CNRL cluster paths
projectDir = '/crnldata/socialhealth/projects/2024_PainGen/';
toolboxdir = fullfile(projectDir, 'Toolboxes');
dataDir = '/crnldata/socialhealth/data/2024_SCEBL_GEN/Imaging_Data/analyses/model3_generalization_ST';
saveDir = '/crnldata/socialhealth/projects/2024_PainGen/results/imaging/STmodel3_cues_vifsQC';

% Add Toolboxes to paths
addpath(genpath(fullfile(toolboxdir, 'RobustToolbox')));
addpath(genpath(fullfile(toolboxdir, 'Neuroimaging_Pattern_Masks')));
addpath(genpath(fullfile(toolboxdir, 'CanlabCore')));
addpath(genpath(fullfile(toolboxdir, 'MediationToolbox')));
addpath(genpath(fullfile(toolboxdir, 'spm12')));

% Initialize variables
subjects = struct();
subjects.vif_thresh = 4;

% Except SPM/external for compatibility with Canlab
% toolbox/mediation 
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

%%  Functionnal images extraction
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
subjects.name = cell(1, numel(subjectsDir));
for p = 1:numel(subjectsDir)
    [~, subname] = fileparts(subjectsDir(p));
    subjects.name{p} = subname;
end 

%subjects.name = subjects.name(1:3); % LOCAL TEST

% Extract specific contrast files from SPM.mat, and saves vifs in subjects struct
% from SPM.mat files, extract contrast of interest, compute vifs/contrast
func_imgs = cell(1, numel(subjects.name)); % Initialize func_imgs as a cell array
subjects.vifs = cell(1, numel(subjects.name));
% get func files in a cell/sub
for sub = 1:numel(subjects.name)
    contrastFiles = {}; % Initialize contrastFiles as an empty cell array
    subFiles = fullfile(dataDir, char(subjects.name{sub}));
    
    SPM = load(fullfile(subFiles, 'SPM.mat')); % Load SPM.mat for the subject
    table = struct2table(SPM.SPM.Vbeta);
    
    cont_names = []
    % Vifs;
    vifRes = scn_spm_design_check(subFiles);
    vifnames = [];
    vifs = [];
    % Loop through each contrast
    for cont = 1:numel(table.descrip)
        % Check if the description of the contrast contains 'allPains_trials0'
        if contains(table.descrip{cont}, 'AllPains_trial0')

            cont_names = [cont_names, table.descrip{cont}]; % get contrast name to find corresponding contrast to look at vifs after
            
            % add the image path to contrastFiles
            add_img = table.fname{cont};
            img_path = fullfile(dataDir, subjects.name{sub}, add_img);
            contrastFiles = [contrastFiles, img_path]; % Append image path to contrastFiles

            %vifs
            vifnames = [vifnames; cellstr(vifRes.name{cont})];
            vifs = [vifs; vifRes.allvifs(cont)];
        end
    end
    
    subjects.vifs{sub} = vifs;
    subjects.cont{sub} = vifnames;
    func_imgs{sub} = char(contrastFiles'); % Assign contrastFiles to the corresponding subject in func_imgs
end

%% look for vifs outlier and save in subject struct
total_outlier = 0;
cont_count = 0;
for subject_idx = 1:numel(subjects.name)
    
    sub_vifs = subjects.vifs{subject_idx};
    is_out = zeros(1,numel(sub_vifs));
    
    for j=1:numel(sub_vifs) % check each contrast vifs and save True at idx j if > thresh
        if sub_vifs(j) > subjects.vif_thresh
            is_out(j) = 1;
            total_outlier = total_outlier + 1; % to know total num of outlier
        end
    end
    cont_count = cont_count + j; % save number of total iterations
    disp(['For subject ', subjects.name{subject_idx}, ', contrasts at indices ', num2str(find(is_out == 1)), ' were above the threshold of ', num2str(subjects.vif_thresh)]);
    subjects.outvifs{subject_idx} = is_out;
end
disp(['Total number of VIF outliers at threshold ', num2str(subjects.vif_thresh), ' is ', num2str(total_outlier),...
    '/', num2str(cont_count), ' total contrasts (',num2str(total_outlier/cont_count*100), '%)' ]);

cd(saveDir)
save ('subjectsVifs.mat', "subjects")

%% plotvifs

figure;
hold on; 
for subject_idx = 1:numel(subjects.name)
    
    sub_vifs = subjects.vifs{subject_idx};

    % Connect the data points with lines
    plot(1:numel(sub_vifs), sub_vifs, 'LineWidth', 1);
    scatter(1:numel(sub_vifs), sub_vifs, 50, 'filled');
end
hold off; % Disable hold after overlaying plots
yline(subjects.vif_thresh, '--', 'Color', 'r', 'LineWidth', 1.5); 

ylim([0 10]);
title(['VIFs per cue ST contrast for all ', num2str(numel(subjects.name)), ' subjects']);
xlabel('Cue-evoked contrasts');
ylabel('VIF value');
annotation('textbox', [0.1, 0.05, 0.8, 0.1], 'String', {['Total number of VIF outliers at threshold = ',...
    num2str(subjects.vif_thresh), ' is ', num2str(total_outlier), '/', num2str(cont_count),' (',num2str(total_outlier/cont_count*100), '%) contrasts images']}, 'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'FontSize', 12);


% %% Qual check
% data = fmri_data(func_imgs)
% 
% plot(data)
% 
% 
% wh_cols = scn_spm_get_events_of_interest(pwd, 'events_only')
% for each sub:
% 
%     out = scn_spm_design_check(subFiles, 'events_only')
% 
%     target_vifs = 
% spm_general_hist(func_imgs, func_imgs, 'all')


