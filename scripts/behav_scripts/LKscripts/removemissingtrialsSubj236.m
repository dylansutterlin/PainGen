%% correct for missing trials of Subject 236
% in order to have same indexes for behavioral ratings and beta images (92
% trials instead of 96, because run 4 was aborted 4 runs before the end
% LK, May 2016

clc
clear all
close all

%% load behavioral data

cd ('/Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Behavioral/Results')
load SCEBLmri_data_FINAL_N36;

fs = fieldnames(SCEBLmri_data);
for f = 7:numel(fs)
    SCEBLmri_data.(fs{f}){36} = SCEBLmri_data.(fs{f}){36}([1:60 65:96]);   
end

%% load NPS expression data

load /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/pain_exp_38.mat
load /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/painvifs.mat

% remove high vif trials from single trial NPS data
for s = 1:38
    SCEBLdata.high_vif_trials_idx{1,s} = find(painvifs{s} > 2.5);
    pain_exp{s}(SCEBLdata.high_vif_trials_idx{1,s}) = NaN;
end

%% load physio data

glmfiles = filenames('/Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Physio/Data/SCEBL_MRI*/Biopac/GLM_SCR_onlypain_singletrials/gsr.mat');

for subj = 1:numel(glmfiles)
    load(glmfiles{subj})
    if subj ==36
        betas_all{subj} = glm.beta(1:92);
    else
        betas_all{subj} = glm.beta(1:96);
    end
 end

%% Take only medium temperature (49) trials to control completely for temperature

for subj = 1:38 
    
    T49.inds{subj} = find(SCEBLmri_data.temp{subj} == 49);
    
    T49.trialinblock{subj} = SCEBLmri_data.blocktrial{subj}(T49.inds{subj}); 
    
    T49.cues{subj} = SCEBLmri_data.cues{subj}(T49.inds{subj}); 
    T49.social{subj} = SCEBLmri_data.social{subj}(T49.inds{subj}); 
        
    T49.cues_trial{subj} = [SCEBLmri_data.cues{subj}(T49.inds{subj}) meancenter(T49.trialinblock{subj})]; 
    T49.social_trial{subj} = [SCEBLmri_data.social{subj}(T49.inds{subj}) meancenter(T49.trialinblock{subj})]; 
       
    T49.exp{subj} = SCEBLmri_data.exp{subj}(T49.inds{subj}); 
    T49.pain{subj} = SCEBLmri_data.pain{subj}(T49.inds{subj}); 
    
    T49.gsr{subj} = betas_all{subj}(T49.inds{subj});
    
    T49.nps{subj} = pain_exp{subj}(T49.inds{subj});
    T49.nps_trialinblock{subj} = [pain_exp{subj}(T49.inds{subj}) T49.trialinblock{subj}];
    
    T49.mfpe{subj} = SCEBLmri_data.mfpe{subj}(T49.inds{subj}); 
                   
end

save T49_236corr T49
save SCEBLmri_data_236corr SCEBLmri_data

