%% Analyze behavioral data SCEBL
% by LK Oct 31, 2013
% This script is based directly on the txt-file (not edat file!!) and
% slightly different than the one used for the pre-pilots

addpath('/Users/Leonie/Documents/Tools/scripts/General scripts')

clc
close all
clear all
cd C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Behavioral/DATA

subjects = filenames('SCEBL_MRI2*');

% cm = [0 0 0; .2 .2 .2; 0 0 .5; 0 0 .9; .5 0 0; 1 0 0; .5 0 .5; 1 0 .8];
cm = [.2 .2 .2; .4 .4 .4; 0 0 .6; 0 0 .9; .7 0 0; 1 0 0; .5 0 .2; 1 0 .4];
cms = [.1 .1 .1; 0 0 .7; .75 0 0; .75 0 .75];
cmtemp = [.4 .1 .1; .6 .3 .3; 1 .4 .4];



for sub = 1:numel(subjects)

    cd('C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Behavioral/DATA/')

    % put current subject folder here:
    subject = (subjects{sub})

    cd (subject)
    cd ('3_Main')
    
   f = filenames('L*.xlsx', 'char');
    
   SCEBLmri_data.subjname{sub} = (subjects{sub});
   SCEBLmri_data.cb{sub} = xlsread(f, [f(1:4), '.txt'], 'D3');
   [DDD SCEBLmri_data.hc{sub}] = xlsread(f, [f(1:4), '.txt'], 'M3');
   [DDD2 SCEBLmri_data.lc{sub}] = xlsread(f, [f(1:4), '.txt'], 'N3'); clear DDD DDD2
   SCEBLmri_data.occid{sub} = xlsread(f, [f(1:4), '.txt'], 'B3');

   SCEBLmri_data.forcedchoice{sub} = xlsread(f, [f(1:4), '.txt'], 'AG3:AG104');
   SCEBLmri_data.forcedchoice{sub} = SCEBLmri_data.forcedchoice{sub}(~isnan(SCEBLmri_data.forcedchoice{sub}));

   dummy = xlsread(f, [f(1:4), '.txt'], 'CE3:CE103'); 
   SCEBLmri_data.blocktrial{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);

   SCEBLmri_data.block{sub} = [repmat(1, 16, 1); repmat(2, 16, 1); repmat(3, 16, 1); repmat(4, 16, 1); repmat(5, 16, 1); repmat(6, 16, 1)];

   dummy = xlsread(f, [f(1:4), '.txt'], 'CF3:CF103'); 
   SCEBLmri_data.temp{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);

   dummy = xlsread(f, [f(1:4), '.txt'], 'CC3:CC103'); 
   SCEBLmri_data.condition{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);

%    dummy = xlsread(f, [f(1:4), '.txt'], 'AJ3:AJ103'); 
%    SCEBLmri_data.condname{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'AT3:AT103'); 
   SCEBLmri_data.RTexp{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'BO3:BO103'); 
   SCEBLmri_data.RTpain{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'AP3:AP103'); 
   SCEBLmri_data.exp{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'BK3:BK103'); 
   SCEBLmri_data.pain{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
   
   SCEBLmri_data.mfpe{sub} = SCEBLmri_data.pain{sub}-SCEBLmri_data.exp{sub};
   
   SCEBLmri_data.exp_mfpe{sub} = [SCEBLmri_data.exp{sub}, SCEBLmri_data.mfpe{sub}];
   
   
    % assign 1 to cue-high trials and -1 to cue-low trials
    cl_trials = find(SCEBLmri_data.condition{sub} < 5);
    ch_trials = find(SCEBLmri_data.condition{sub} > 4);
    ctrials(cl_trials) = -1;
    ctrials(ch_trials) = 1;
    SCEBLmri_data.cues{sub} = ctrials';

    % assign 1 to social-high trials and -1 to social-low trials  
    sl_trials = find(SCEBLmri_data.condition{sub} == 1 | SCEBLmri_data.condition{sub} == 2 | SCEBLmri_data.condition{sub} == 5 | SCEBLmri_data.condition{sub} == 6);
    sh_trials = find(SCEBLmri_data.condition{sub} == 3 | SCEBLmri_data.condition{sub} == 4 | SCEBLmri_data.condition{sub} == 7 | SCEBLmri_data.condition{sub} == 8);
    strials(sl_trials) = -1;
    strials(sh_trials) = 1;
    SCEBLmri_data.social{sub} = strials';


    % assign 1 to conflict trials and -1 to noconflict trials  
    SCEBLmri_data.conflict{sub} = -1 * ones(96,1);
    SCEBLmri_data.conflict{sub}(find(SCEBLmri_data.condition{sub} > 2 & SCEBLmri_data.condition{sub} < 7)) = 1;
    
    SCEBLmri_data.cues_soc_int{sub} = [SCEBLmri_data.cues{sub}, SCEBLmri_data.social{sub}, SCEBLmri_data.conflict{sub}];

    % exact social stim values
  
    clear socmean socstd
    
    [dummy2 dummy] = xlsread(f, [f(1:4), '.txt'], 'BZ3:BZ103'); clear dummy2
    socstims = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);

%     for t = 1:96
%         stimname = char(socstims(t));
%         socmean(t) = str2num(stimname(end-13:end-11))/1000;
%         socstd(t) = str2num(stimname(end-6:end-4))/1000;
%     end
% 
%     SCEBLmri_data.socmeans{sub} = socmean;
%     SCEBLmri_data.socSTDs{sub} = socstd;
%     SCEBLmri_data.MSTDint{sub} = meancenter(socmean) .* meancenter(socstd);
    

    %% Plot by temperature

    figure('Color', [1 1 1], 'Name', (subjects{sub}));
    
    subplot(3,3,1)
    scatter(SCEBLmri_data.temp{sub}, SCEBLmri_data.pain{sub});
    xlim([47.5 50.5])
        
    subplot(3,3,4)
    plot(find(SCEBLmri_data.temp{sub} == 48), SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 48), 'Color', [.4 .1 .1], 'LineWidth', 1.5); hold on;
    plot(find(SCEBLmri_data.temp{sub} == 49), SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 49), 'Color', [.6 .3 .3], 'LineWidth', 1.5); hold on;
    plot(find(SCEBLmri_data.temp{sub} == 50), SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 50), 'Color', [1 .4 .4],  'LineWidth', 1.5); hold on;
    legend('48C', '49C', '50C');

    subplot(3,3,7)
    bar([48,49,50], [nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 48)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 49)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 50))]);
    xlim([47.5 50.5])

    %% Plot by condition
    
    clsltrials = find(SCEBLmri_data.condition{sub} == 1 | SCEBLmri_data.condition{sub} == 2);
    clshtrials = find(SCEBLmri_data.condition{sub} == 3 | SCEBLmri_data.condition{sub} == 4);
    chsltrials = find(SCEBLmri_data.condition{sub} == 5 | SCEBLmri_data.condition{sub} == 6);
    chshtrials = find(SCEBLmri_data.condition{sub} == 7 | SCEBLmri_data.condition{sub} == 8);

    %% EXPT Timecourse

%     figure('Color',[1 1 1],  'Name', (subjects{sub}));
    
    subplot(2,3,2)
    plot(clsltrials, SCEBLmri_data.exp{sub}(clsltrials), 'Color', cms(1,:), 'LineWidth',1.5); hold on; 
    plot(clshtrials, SCEBLmri_data.exp{sub}(clshtrials), '--', 'Color', cms(2,:), 'LineWidth',1.5); hold on; 
    plot(chsltrials, SCEBLmri_data.exp{sub}(chsltrials), '--', 'Color', cms(3,:), 'LineWidth',1.5); hold on; 
    plot(chshtrials, SCEBLmri_data.exp{sub}(chshtrials), 'Color', cms(4,:), 'LineWidth',1.5); hold on; 
    legend('CLSL', 'CLSH', 'CHSL', 'CHSH');
    title('Expectation ratings ALL trials','FontSize',16, 'FontName','Arial');
    
    % Barplots Expectations
    
    subplot(2,3,5)
    if sub == 34
        barplot_colored([SCEBLmri_data.exp{sub}(clsltrials), SCEBLmri_data.exp{sub}(clshtrials), [SCEBLmri_data.exp{sub}(chsltrials); NaN; NaN], [SCEBLmri_data.exp{sub}(chshtrials); NaN]])
    else 
        barplot_colored([SCEBLmri_data.exp{sub}(clsltrials), SCEBLmri_data.exp{sub}(clshtrials), SCEBLmri_data.exp{sub}(chsltrials), SCEBLmri_data.exp{sub}(chshtrials)])
    end
    set(gca,'XTickLabel',{'CLSL', 'CLSH', 'CHSL', 'CHSH'});
    title('Expectation Ratings ALL trials','FontSize',16, 'FontName','Arial');
    colormap(cms)

    
    %% timecourse medium temp PAIN ratings by condition

%     figure('Color',[1 1 1],  'Name', (subjects{sub}));
    
    subplot(3,3,3)
    plot(clsltrials, SCEBLmri_data.pain{sub}(clsltrials), 'Color', cms(1,:), 'LineWidth',1.5); hold on; 
    plot(clshtrials, SCEBLmri_data.pain{sub}(clshtrials), '--', 'Color', cms(2,:), 'LineWidth',1.5); hold on; 
    plot(chsltrials, SCEBLmri_data.pain{sub}(chsltrials), '--', 'Color', cms(3,:), 'LineWidth',1.5); hold on; 
    plot(chshtrials, SCEBLmri_data.pain{sub}(chshtrials), 'Color', cms(4,:), 'LineWidth',1.5); hold on; 
    legend('CLSL', 'CLSH', 'CHSL', 'CHSH');
    title('PAIN ratings all heat trials','FontSize',16, 'FontName','Arial');

    subplot(3,3,6)
    if sub == 34 
        barplot_colored([SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 2), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 4), [SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 5); NaN], SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 7)])
    else
        barplot_colored([SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 2), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 4), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 5), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 7)])
    end
    set(gca,'XTickLabel',{'CLSL49', 'CLSH49', 'CHSL49', 'CHSH49'});
    title('PAIN Ratings medium heat','FontSize',16, 'FontName','Arial');
    colormap(cms)

    subplot(3,3,9)
    if sub == 34
        barplot_colored([SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 1), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 3), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 2), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 4), [SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 5);NaN], SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 7), [SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 6);NaN], [SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 8);NaN]])
    else
        barplot_colored([SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 1), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 3), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 2), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 4), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 5), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 7), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 6), SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 8)])
    end
    set(gca,'XTickLabel',{'CLSL48', 'CLSL49', 'CLSH48', 'CLSH49', 'CHSL49', 'CHSL50', 'CHSH49', 'CHSH50'});
    title('PAIN Ratings all','FontSize',16, 'FontName','Arial');
    colormap(cm)

    %%
    painmeans_all(sub,:) = [nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 1)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 3)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 2)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 4)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 5)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 7)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 6)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 8))];
    expemeans_all(sub,:) = [nanmean(SCEBLmri_data.exp{sub}(clsltrials)), nanmean(SCEBLmri_data.exp{sub}(clshtrials)), nanmean(SCEBLmri_data.exp{sub}(chsltrials)), nanmean(SCEBLmri_data.exp{sub}(chshtrials))];
    pain_all_temp(sub,:) = [nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 48)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 49)), nanmean(SCEBLmri_data.pain{sub}(SCEBLmri_data.temp{sub} == 50))];
    painRT_all(sub,:) = [nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 1)), nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 3)), nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 2)), nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 4)), nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 5)), nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 7)), nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 6)), nanmean(SCEBLmri_data.RTpain{sub}(SCEBLmri_data.condition{sub} == 8))];
    expeRT_all(sub,:) = [nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 1)), nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 3)), nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 2)), nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 4)), nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 5)), nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 7)), nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 6)), nanmean(SCEBLmri_data.RTexp{sub}(SCEBLmri_data.condition{sub} == 8))];

    if sub == 34
        exp_tc.clsl(sub,:) = SCEBLmri_data.exp{sub}(clsltrials);
        exp_tc.clsh(sub,:) = SCEBLmri_data.exp{sub}(clshtrials);
        exp_tc.chsl(sub,:) = [SCEBLmri_data.exp{sub}(chsltrials); NaN; NaN];
        exp_tc.chsh(sub,:) = [SCEBLmri_data.exp{sub}(chshtrials); NaN];
        pain49_tc.clsl(sub,:) = SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 2);
        pain49_tc.clsh(sub,:) = SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 4);
        pain49_tc.chsl(sub,:) = [SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 5); NaN];
        pain49_tc.chsh(sub,:) = SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 7);

    else
        exp_tc.clsl(sub,:) = SCEBLmri_data.exp{sub}(clsltrials);
        exp_tc.clsh(sub,:) = SCEBLmri_data.exp{sub}(clshtrials);
        exp_tc.chsl(sub,:) = SCEBLmri_data.exp{sub}(chsltrials);
        exp_tc.chsh(sub,:) = SCEBLmri_data.exp{sub}(chshtrials);
        pain49_tc.clsl(sub,:) = SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 2);
        pain49_tc.clsh(sub,:) = SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 4);
        pain49_tc.chsl(sub,:) = SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 5);
        pain49_tc.chsh(sub,:) = SCEBLmri_data.pain{sub}(SCEBLmri_data.condition{sub} == 7);

    end
    
    
  % input('Press any key to continue to next subject');
%     if sub < 30
%         close all
%     end
    
end

% input('Press any key to continue to grand average')


%% Temp effects

figure('Color', [1 1 1])
barplot_colored(pain_all_temp, 'within'); hold on
plot([1 2 3], pain_all_temp, 'ko'); hold on
set(gca,'XTickLabel',{'48C', '49C', '50C'});
title('Pain ratings by temperatures','FontSize',16, 'FontName','Arial');
colormap(cmtemp)

%% add here Todd's pain ratings in red

% plot([1 2 3], [20 40 60], 'r*'); hold on

% output mean and 95% CI values for temperatures
temp_means = mean(pain_all_temp)
temp_ste = ste(pain_all_temp)
temp_95ci_l = temp_means - 2.001 * temp_ste
temp_95ci_h = temp_means + 2.001 * temp_ste


%%
figure('Color', [1 1 1]) 
subplot(2,1,1)
barplot_colored(expemeans_all, 'within')
set(gca,'XTickLabel',{'CLSL', 'CLSH', 'CHSL', 'CHSH'});
title('Expectation Ratings','FontSize',16, 'FontName','Arial');
colormap(cms)

subplot(2,1,2)
barplot_colored(painmeans_all(:,3:6), 'within')
set(gca,'XTickLabel',{'T49 CLSL', 'T49 CLSH', 'T49 CHSL', 'T49 CHSH'});
title('Pain Ratings medium temperatures','FontSize',16, 'FontName','Arial');
colormap(cms)

figure('Color', [1 1 1]) 
barplot_colored(painmeans_all, 'within')
set(gca,'XTickLabel',{'T48 CLSL', 'T48 CLSH', 'T49 CLSL', 'T49 CLSH', 'T49 CHSL', 'T49 CHSH', 'T50 CHSL', 'T50 CHSH'});
title('Pain Ratings ALL trials and temperatures','FontSize',16, 'FontName','Arial');
colormap(cm)

figure('Color', [1 1 1]) 
subplot(2,1,1)
barplot_colored(expeRT_all, 'within')
set(gca,'XTickLabel',{'T48 CLSL', 'T48 CLSH', 'T49 CLSL', 'T49 CLSH', 'T49 CHSL', 'T49 CHSH', 'T50 CHSL', 'T50 CHSH'});
title('Expectation Ratings RTs ALL trials and temperatures','FontSize',16, 'FontName','Arial');
colormap(cm)

subplot(2,1,2)
barplot_colored(painRT_all, 'within')
set(gca,'XTickLabel',{'T48 CLSL', 'T48 CLSH', 'T49 CLSL', 'T49 CLSH', 'T49 CHSL', 'T49 CHSH', 'T50 CHSL', 'T50 CHSH'});
title('Pain Ratings RTs ALL trials and temperatures','FontSize',16, 'FontName','Arial');
colormap(cm)

%% social influence effect for all temperatures

figure('Color', [1 1 1]) 
subplot(1,2,1)
barplot_colored([painmeans_all(:,1), mean([painmeans_all(:,3), painmeans_all(:,5)], 2), painmeans_all(:,7)], 'within')
colormap(cm)
ylim([0 50])

subplot(1,2,2)
barplot_colored([painmeans_all(:,2), mean([painmeans_all(:,4), painmeans_all(:,6)], 2), painmeans_all(:,8)], 'within')
colormap(cm)
ylim([0 50])


%% Time course learning on average

% One figure for time course of expectation and pain ratings

figure('Color', [1 1 1]) ;
subplot(2,1,1)
plot(1:24, nanmean(exp_tc.clsl), 'Color', cm(3,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
plot(1:24, nanmean(exp_tc.clsh), 'Color', cm(4,:), 'LineWidth', 1.5); hold on
plot(1:24, nanmean(exp_tc.chsl), 'Color', cm(5,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
plot(1:24, nanmean(exp_tc.chsh), 'Color', cm(6,:), 'LineWidth', 1.5); hold on
title('Expectation ratings by time', 'FontSize', 16, 'FontName', 'Arial');
xlim([0 25])
legend('clsl', 'clsh', 'chsl', 'chsh');

subplot(2,1,2)
plot(1:12, nanmean(pain49_tc.clsl), 'Color', cm(3,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
plot(1:12, nanmean(pain49_tc.clsh), 'Color', cm(4,:), 'LineWidth', 1.5); hold on
plot(1:12, nanmean(pain49_tc.chsl), 'Color', cm(5,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
plot(1:12, nanmean(pain49_tc.chsh), 'Color', cm(6,:), 'LineWidth', 1.5); hold on
title('Pain ratings medium temp by time', 'FontSize', 16, 'FontName', 'Arial');
xlim([0.5 12.5])
legend('clsl', 'clsh', 'chsl', 'chsh');



%% regression lines for relationship between expectation and pain ratings

figure('Color', [1 1 1])
h = line_plot_multisubject(SCEBLmri_data.exp, SCEBLmri_data.pain); title('Multi-subject scatterplot: expectation ratings - pain ratings','FontSize',16, 'FontName', 'Arial', 'FontWeight', 'Bold')
xlabel('Expectation', 'FontSize', 14, 'FontName','Arial');
ylabel('Pain', 'FontSize', 14, 'FontName','Arial');


% cd C:\Users\Leonie\Dropbox\work\projects\6_SCEBL1_fMRI\Behavioral/Results
% save SCEBLmri_data SCEBLmri_data



%% forced choice analysis

SCEBLmri_data.forcedchoice{34} = 6*ones(6,1);

for sub = 1:numel(subjects)
    forcedchoice_all(sub,:) = SCEBLmri_data.forcedchoice{sub}';
end

forcedchoice_all_bin = forcedchoice_all > 1;

figure('Color', [1 1 1], 'Name', 'Forced Choice time course runs 1-6')
barplot_colored(forcedchoice_all_bin, 'within'); 
title('Forced choice time course')
xlabel('Runs 1-6')
ylabel('Correct choice = 1')

    
% %% for (Tor's) single trial analyses
% %  
% load /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Imaging/Analyses/painvifs.mat
% 
% for s = 1:38
%     SCEBLdata_singletrial.high_vif_trials_idx{1,s} = find(painvifs{s} > 2.5);
% end
% 
% SCEBLdata_singletrial.temp = SCEBLmri_data.temp(1:38);
% SCEBLdata_singletrial.ratings = SCEBLmri_data.pain(1:38);
% SCEBLdata_singletrial.high_vif_descript = 'vif > 2.5';
% SCEBLdata_singletrial.subjects = subjects(1:38)';
% SCEBLdata_singletrial.dat_obj = filenames('painbetas*.mat')';
% 
% 
% % save ('SCEBLdata_singletrial_N38.mat', 'SCEBLdata_singletrial')
%     
% %% for canlab dataset
% 
% SCEBL_fMRI = canlab_dataset
% 
% 
% SCEBL_fMRI.Event_Level.names  = fieldnames(SCEBLmri_data)';
% SCEBL_fMRI.Event_Level.names  = SCEBL_fMRI.Event_Level.names([7:15 17:19])
% 
% for s = 1:38
%     
%     SCEBL_fMRI.Event_Level.data{s} = [SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{1}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{2}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{3}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{4}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{5}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{6}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{7}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{8}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{9}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{10}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{11}){s}, ...
%                                     SCEBLmri_data.(SCEBL_fMRI.Event_Level.names{12}){s}];
% end
% 
% SCEBL_fMRI.Description.Experiment_Name = 'SCEBL_fMRI';
% 
% SCEBL_fMRI.Subj_Level.id = SCEBLmri_data.subjname';
% 
% 
% 
% 
