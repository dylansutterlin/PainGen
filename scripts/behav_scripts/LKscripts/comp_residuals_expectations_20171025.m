%% add fields and do glmfit multilevel to compute residuals for expectations
% LK Oct 2017

clear all
close all
clc

cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Behavioral/Results
load('SCEBLmri_data_236corr.mat')

for s = 1:38
   SCEBLmri_data.cues_soc{s} = [SCEBLmri_data.cues{s}, SCEBLmri_data.social{s}]; 
end

stats = glmfit_multilevel(SCEBLmri_data.exp, SCEBLmri_data.cues_soc, []);


%% compute residual expectations (covarying out social and cue effects)

figure

for s = 1:38
    
    SCEBLmri_data.expfit{s} = stats.first_level.beta(1,s) + stats.first_level.beta(2,s)*SCEBLmri_data.cues{s} + stats.first_level.beta(3,s)*SCEBLmri_data.social{s};
%     subplot(6,7,s); scatter(SCEBLmri_data.expfit{s}, SCEBLmri_data.exp{s});
    
    SCEBLmri_data.expresid{s} = SCEBLmri_data.exp{s}-SCEBLmri_data.expfit{s};
    subplot(6,7,s); plot(SCEBLmri_data.expresid{s})
    
end
    
%% for medium temp trials

load('T49_236corr.mat')

for s = 1:38
    
    T49.expresid{s} = SCEBLmri_data.expresid{s}(T49.inds{s});
    
end