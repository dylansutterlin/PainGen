

% script to extract output from LKscript/analyze_* main, simrate and gen
% Load its outputs/var to perform main effects
% Save struct as dataframe in ~/results/ for further stats in indendent
% scripts


clc
close all
clear all
projectDir = '/home/dsutterlin/projects/genPain/';
dataDir = fullfile(projectDir, 'DATA/Behavioral/')
saveDir = fullfile(projectDir, 'results/behavioral')

cd (saveDir);
load SCEBLmri_learndata_FINAL_N36.mat
load simmat.mat
%load simrate.mat
load SCEBLmri_data_gen36.mat 

smsim = xlsread('smDistance-wingfield2022.xlsx') % sensorimotor cosine sim
% stim order : i=[dog horse cow car train truck]

            
cd (dataDir);
subjects = filenames('SCEBL_MRI2*');

%%%%%%%%%%%%%%%
% Learning task
fieldNames = fieldnames(SCEBLmri_data);

% Specify the fields you want to process
targetFields = {'blocktrial', 'block', 'temp', 'condition', 'RTexp', 'RTpain', 'exp', 'pain', 'mfpe', 'cues', 'social', 'conflict'};

%% Strict categ (1/-1) theoretical similarity matrix
A = [1 -1 -1 -1 -1 -1;
     -1 1 -1 -1 -1 -1;
     -1 -1 1 -1 -1 -1;
     -1 -1 -1 1 -1 -1;
     -1 -1 -1 -1 1 -1;
     -1 -1 -1 -1 -1 1];
categSimmat = repmat(A, 1, 3) % 6 (learning stim) x 18 (gen stim)

figure('Color', [1 1 1], 'Name','Strict conceptual similarity ratings')
imagesc(categSimmat);
ylabel('Learning cues');
yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);

xIndexNames = {'GCA1'; 'GCA2'; 'GCA3'; 'GCV1'; 'GCV2'; 'GCV3';...
           'GPA1'; 'GPA2'; 'GPA3'; 'GPV1'; 'GPV2'; 'GPV3'; ...
           'GWA1'; 'GWA2'; 'GWA3'; 'GWV1'; 'GWV2'; 'GWV3'};
xlabel('Generalization cues');
set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);
colorbar

%% Soft categ (1/-1) theoretical similarity matrix
A = [1 1 1 -1 -1 -1;
     1 1 1 -1 -1 -1;
     1 1 1 -1 -1 -1;
     -1 -1 -1 1 1 1;
     -1 -1 -1 1 1 1;
     -1 -1 -1 1 1 1];
softcategSimmat = repmat(A, 1, 3) % 6 (learning stim) x 18 (gen stim)

figure('Color', [1 1 1], 'Name','Soft categ conceptual similarity ratings')
imagesc(softcategSimmat);
ylabel('Learning cues');
yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);
xlabel('Generalization cues');
set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);
colorbar

%% mean model between soft and strict categ (1/-1) theoretical similarity matrix

mean_model = (softcategSimmat + categSimmat) /2
figure('Color', [1 1 1], 'Name','Mean model (two categ) similarity ratings')
imagesc(mean_model);
ylabel('Learning cues');
yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);
xlabel('Generalization cues');
set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);
colorbar
%% Semantic(word2vec; Boulder U) similarity matrix 
A = [
    1       0.483   0.402   0.31    0.216   0.247;
    0.483   1       0.494   0.249   0.25    0.226;
    0.402   0.494   1       0.176   0.174   0.232;
    0.31    0.249   0.176   1       0.34    0.674;
    0.216   0.25    0.174   0.34    1       0.38;
    0.247   0.226   0.232   0.674   0.38    1
];
boulderSimmat = repmat(A, 1, 3) % 6 (learning stim) x 18 (gen stim)

figure('Color', [1 1 1], 'Name','Semantic word2vec (Boulder) conceptual similarity')
imagesc(boulderSimmat);
ylabel('Learning cues');
yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);
xlabel('Generalization cues');
set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);
colorbar

%% Semantic Boulder (LSA) similarity matrix
A = [1      0.483  0.402  0.31   0.216  0.247;
     0.483  1      0.494  0.249  0.25   0.226;
     0.402  0.494  1      0.176  0.174  0.232;
     0.31   0.249  0.176  1      0.34   0.674;
     0.216  0.25   0.174  0.34   1      0.38;
     0.247  0.226  0.232  0.674  0.38   1];
boulderLSASimmat = repmat(A, 1, 3) % 6 (learning stim) x 18 (gen stim)

figure('Color', [1 1 1], 'Name','Semantic-LSA (Boulder) conceptual similarity')
imagesc(boulderLSASimmat);
ylabel('Learning cues');
yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);
xlabel('Generalization cues');
set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);
colorbar
%% SensoriMotor (Landcaster norm) similarity matrix

smSimmat = repmat(smsim, 1, 3) % 6 (learning stim) x 18 (gen stim)
smSimmat = 1-smSimmat;
figure('Color', [1 1 1], 'Name','SensoriMotor conceptual similarity ratings')
imagesc(smSimmat);
ylabel('Learning cues');
yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);
xlabel('Generalization cues');
set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);
colorbar

%% Empirical similarity matrix
lc_sim = squeeze(simmat3(1, :, :))'; % Extract first row and transpose
hc_sim = squeeze(simmat3(2, :, :))'; % Extract second row and transpose

simlearngen = struct('dog', [], 'horse', [], 'cow', [], 'car', [], 'train', [], 'truck', []);
%lowsim = struct('dog', [], 'horse', [], 'cow', [], 'car', [], 'train', [], 'truck', []);

% building two struct for ratings grouped on learning cues
for sub = 1:numel(subjects)
    % get learning cues from this sub
    lc = SCEBLmri_gendata.lowcue{sub};
    hc = SCEBLmri_gendata.highcue{sub};
    
    % ratings for this sub
    rowlc = lc_sim(sub, :);
    rowhc = hc_sim(sub, :);

    % Handling High Conditions
    if contains(hc, 'dog')
        simlearngen.dog = [simlearngen.dog; rowhc];
    elseif contains(hc, 'horse')
        simlearngen.horse = [simlearngen.horse; rowhc];
    elseif contains(hc, 'cow')
        simlearngen.cow = [simlearngen.cow; rowhc];
    elseif contains(hc, 'car')
        simlearngen.car = [simlearngen.car; rowhc];
    elseif contains(hc, 'train')
        simlearngen.train = [simlearngen.train; rowhc];
    elseif contains(hc, 'truck')
        simlearngen.truck = [simlearngen.truck; rowhc];
    end

end

% Computing similarity matrix
empiricalSimmat = nan(6, 18);
learncue = fieldnames(simlearngen);
for i = 1:numel(learncue)
    fieldName = learncue{i};
    medianVector = nanmedian(simlearngen.(fieldName));
    % Place this vector at the corresponding row in empiricalsimmat
    empiricalSimmat(i, :) = medianVector;
end
    
figure('Color', [1 1 1])
imagesc(empiricalSimmat);
title('Empirical similarity ratings')
ylabel('Learning cues');
yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);
xlabel('Generalization cues');
set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);  
colorbar

%% Ponderate empirical similarity matrix with cosine sim from sensori-m distance
% https://embodiedcognitionlab2.shinyapps.io/sensorimotordistance/

smsimmat = nan(6, 18);
learncue = fieldnames(simlearngen);
for i = 1:numel(learncue);
    fieldName = learncue{i};
    
    medianVector = nanmedian(simlearngen.(fieldName));
    medianVector = medianVector/100; % normal 0-1
    % Extract columns i from sensorim mat 3 times (i=1 dog-dog, dog-horse,
    % dog-cow, dog-car, dog-train, dog-truck.. X3 to match gen cues
    % modalities
    sm = smsim(:, i);
    sm = repmat(sm, 1, 3);
    sm = sm(:)'; %flatten and transpose
    sm = 1-sm; % 1 is now more similar than 0
    % Multiply element-wise with the median vector
    wmedianVector = medianVector .* sm;
    wmedianVector = wmedianVector .* 100;
    %wmedianVector = 1-medianVector
    % Place the resulting vector at the corresponding row in smsimmat
    smsimmat(i, :) = wmedianVector;
end

% figure('Color', [1 1 1], 'Name','Weighted (sensoriM) Empirical conceptual similarity ratings')
% imagesc(smsimmat);
% ylabel('Learning cues');
% yIndexNames = {'dog', 'horse', 'cow', 'car', 'train', 'truck'}
% set(gca, 'YTick', 1:numel(yIndexNames), 'YTickLabel', yIndexNames);
% xlabel('Generalization cues');
% set(gca, 'XTick', 1:numel(xIndexNames), 'XTickLabel', xIndexNames);
% colormap(flipud(parula));  
% colorbar

%% Learning task : Convert fields from learn struct to a . csv (Sub X field) df
df_learnData = []
for sub = 1:numel(subjects);
    
    concatenatedData = [];

    for i = 1:numel(targetFields);
        
        fieldName = targetFields{i};
        fieldContent = SCEBLmri_data.(fieldName);
        cellData = cell2mat(fieldContent(sub));

        % Concatenate the data horizontally (sub 1 ../ sub 2../)
        concatenatedData = [concatenatedData, cellData];
        size(concatenatedData)
    end
    % Store the concatenated data for this subject
    df_learnData = [df_learnData; concatenatedData]; 
end

% Vertically concatenate the data from all fields, i.e. transform-->df
finalgenData = array2table(df_learnData);
finalgenData.Properties.VariableNames = targetFields

cd (saveDir)
writetable(finalgenData, 'SCEBLmri_Learndata_TxT_N36.csv')

% Display the size of the final data array
%disp(['Size of the final data array: ', mat2str(size(finalData))]);

%%%%%%%%%%%%%%%
%% GEN data to dataframe + compute similarity ratings for each trial

fieldNames = fieldnames(SCEBLmri_data);
allData = []

for sub = 1:numel(subjects)
    % Loop to create a new variable in SCEBLmri_gendata, for each trials,
    % that contains theoretical sim ratings (betweem learning cue and gen
    % cue)

    % SIMILARITY RATING PER TRIAL for each specified simMat
    % ('theotrialSim','emptrialSim', 'smweightTrialSim'...)
    categSimmat_rt = cell(1,0)
    softcategSimmat_rt = cell(1,0);
    mean_model_rt = cell(1,0);
    smSimmat_rt = cell(1,0);
    boulderLSASimmat_rt = cell(1,0)
    boulderSimmat_rt = cell(1,0)
    empiricalSimmat_rt = cell(1,0);
   
    
    % Algo to assign similarity value for each trial based on different
    % simMats : Takes CS+ learning sim - CS- learning sim to compute trial
    % i sim score
    for i = 1:numel(SCEBLmri_gendata.condition{sub});
        condj = SCEBLmri_gendata.condition{sub}(i);
        lc = SCEBLmri_gendata.lowcue{sub};
        hc = SCEBLmri_gendata.highcue{sub}; % img name used for learning
        csTrial = SCEBLmri_gendata.cscat{sub}(i); % -1 or 1
        %truecs = SCEBLmri_gendata.AVcat{sub}(i);

            % low cue index
            if contains(lc, 'dog')
                lowi = 1;
            elseif contains(lc, 'horse')
                lowi = 2;
            elseif contains(lc, 'cow')
                lowi = 3;
            elseif contains(lc, 'car')
                lowi = 4;
            elseif contains(lc, 'train')
                lowi = 5;
            elseif contains(lc, 'truck')
                lowi = 6;
            end
            % High cue index
            if contains(hc, 'dog')
                highi = 1;
            elseif contains(hc, 'horse')
                highi = 2;
            elseif contains(hc, 'cow')
                highi = 3;
            elseif contains(hc, 'car')
                highi = 4;
            elseif contains(hc, 'train')
                highi = 5;
            elseif contains(hc, 'truck')
                highi = 6;
            end

        %Extract sim rating fr this trial for each similarity models
        %(stored in simMats)
        
        diff_smr = categSimmat(highi,condj) - categSimmat(lowi,condj)
        categSimmat_rt = [categSimmat_rt, diff_smr];

        diff_smr = categSimmat(highi,condj) - categSimmat(lowi,condj)
        softcategSimmat_rt = [softcategSimmat_rt, diff_smr];

        diff_smr = mean_model(highi,condj) - mean_model(lowi,condj)
        mean_model_rt = [mean_model_rt, diff_smr];

        diff_smr = smSimmat(highi,condj) - smSimmat(lowi,condj)
        smSimmat_rt = [smSimmat_rt, diff_smr];

        diff_smr = boulderLSASimmat(highi,condj) - boulderLSASimmat(lowi,condj)
        boulderLSASimmat_rt = [boulderLSASimmat_rt, diff_smr];
        
        diff_smr = boulderSimmat(highi,condj) - boulderSimmat(lowi,condj)
        boulderSimmat_rt = [boulderSimmat_rt, diff_smr];
        
        diff_smr = empiricalSimmat(highi,condj) - empiricalSimmat(lowi,condj)
        empiricalSimmat_rt = [empiricalSimmat_rt, diff_smr];
        
        
    end
    
    SCEBLmri_gendata.categSimr{sub} = cell2mat(categSimmat_rt)';
    SCEBLmri_gendata.softcategSimr{sub} = cell2mat(softcategSimmat_rt)';
    SCEBLmri_gendata.mean_modelSimr{sub} = cell2mat(mean_model_rt)';
    SCEBLmri_gendata.smSimr{sub} = cell2mat(smSimmat_rt)';
    SCEBLmri_gendata.boulderLSASimr{sub} = cell2mat(boulderLSASimmat_rt)' ;
    SCEBLmri_gendata.boulderSimr{sub} = cell2mat(boulderSimmat_rt)';
    SCEBLmri_gendata.empiricalSimr{sub} = cell2mat(empiricalSimmat_rt)';
end


% Specify fields that have sub x trials format + simratings that are
% created in loop under
targetFields = {'sub_TxTgrp', 'blocktrial', 'condition', 'AVcat',...
    'modality', 'cscat', 'sim_lc', 'sim_hc', 'sim_diff', 'RTpain',...
    'pain','empiricalSimr', 'categSimr', 'softcategSimr', 'mean_modelSimr', 'smSimr',...
    'boulderLSASimr', 'boulderSimr'};
  
for sub = 1:numel(subjects)
    % Initialize an empty matrix to store the concatenate   d data for this subject
    concatenatedData = [];
    % Iterate through each field
    for ii = 1:numel(targetFields);
        
        fieldName = targetFields{ii};

        % Access the content of the field (1x36 cell)
        fieldContent = SCEBLmri_gendata.(fieldName);
        cellData = cell2mat(fieldContent(sub));

        % Concatenate the data horizontally
        concatenatedData = [concatenatedData, cellData];
        size(concatenatedData);
    end
    
    % Store the concatenated data for this subject
    allData = [allData; concatenatedData];
    
end

% Vertically concatenate the data from all fields
finalData = array2table(allData);
finalData.Properties.VariableNames = targetFields

cd (saveDir)
writetable(finalData, 'SCEBLmri_Gendata_TxT_N36.csv')
save('finalPreproc_simGen_data.mat', 'SCEBLmri_gendata')


%% Correl between similarity ratings
colnames = {'empiricalSimr', 'categSimr', 'softcategSimr', 'mean_modelSimr', 'smSimr',...
    'boulderLSASimr', 'boulderSimr'}
plot_correlation_matrix(finalData(:,colnames), 'names', colnames);
title('Diff (CS+ - CS-) simRatings/trials correlation for simModels');
drawnow, snapnow

figure; plotmatrix(table2array(finalData(:,colnames)));
drawnow, snapnow

colnames = {'blocktrial', 'condition','cscat', 'pain','empiricalSimr', 'mean_modelSimr', 'boulderLSASimr'}
plot_correlation_matrix(finalData(:,colnames), 'names', colnames);
drawnow, snapnow


%% Correlation between empirical similarity mat and other sim models
% Sim matrices to correl with empiricalsimmat
simMats = {'categSimmat', categSimmat; 'softcategSimmat', softcategSimmat; 'mean_model', mean_model; 'smSimmat', smSimmat; 'boulderLSASimmat', boulderLSASimmat; 'boulderSimmat', boulderSimmat};
correlationTable = table();

for i = 1:size(simMats, 1)
    matrixName = simMats{i, 1};
    matrixData = simMats{i, 2};
    spearmanCorrelation = corr(matrixData(:), empiricalSimmat(:), 'Type', 'Spearman');
    correlationTable.(matrixName) = spearmanCorrelation;
end

% Display the correlation table
disp(correlationTable);

