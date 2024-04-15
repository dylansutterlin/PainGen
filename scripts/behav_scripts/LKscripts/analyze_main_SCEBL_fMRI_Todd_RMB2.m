%% Analyze behavioral data SCEBL
% by LK Oct 31, 2013
% This script is based directly on the txt-file (not edat file!!) and
% slightly different than the one used for the pre-pilots

addpath('/Users/Leonie/Documents/Tools/scripts/General scripts')

clc
close all
cd /Users/Leonie/Documents/Buerocratie/Media/RMB2_Sept2014/

subjects = filenames('*Todd');

cm = [0 0 0; .2 .2 .2; 0 0 .5; 0 0 .9; .5 0 0; 1 0 0; .5 0 .5; 1 0 .8];
cms = [.1 .1 .1; 0 0 .7; .75 0 0; .75 0 .75];
cmtemp = [.4 .1 .1; .6 .3 .3; 1 .4 .4];



sub = 1

    % put current subject folder here:
    subject = (subjects{sub})

    cd (subject)
    cd ('3_Main')
    
    f = filenames('L*.xlsx', 'char')
    
   Todd_data.subjname{sub} = (subjects{sub});
   Todd_data.cb{sub} = xlsread(f, [f(1:4), '.txt'], 'D3');
   Todd_data.occid{sub} = xlsread(f, [f(1:4), '.txt'], 'B3');
   
   Todd_data.forcedchoice{sub} = nonnans(xlsread(f, [f(1:4), '.txt'], 'AG3:AG104'));
    
   dummy = xlsread(f, [f(1:4), '.txt'], 'CE3:CE103'); 
   Todd_data.blocktrial{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
   
   dummy = xlsread(f, [f(1:4), '.txt'], 'CF3:CF103'); 
   Todd_data.temp{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
   
   dummy = xlsread(f, [f(1:4), '.txt'], 'CC3:CC103'); 
   Todd_data.condition{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
%    dummy = xlsread(f, [f(1:4), '.txt'], 'AJ3:AJ103'); 
%    SCEBLmri_data.condname{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'AT3:AT103'); 
   Todd_data.RTexp{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'BO3:BO103'); 
   Todd_data.RTpain{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'AP3:AP103'); 
   Todd_data.exp{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
 
   dummy = xlsread(f, [f(1:4), '.txt'], 'BK3:BK103'); 
   Todd_data.pain{sub} = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
   
    % assign 1 to cue-high trials and -1 to cue-low trials
    cl_trials = find(Todd_data.condition{sub} < 5);
    ch_trials = find(Todd_data.condition{sub} > 4);
    ctrials(cl_trials) = -1;
    ctrials(ch_trials) = 1;
   Todd_data.cues{sub} = ctrials';

    % assign 1 to social-high trials and -1 to social-low trials  
    sl_trials = find(Todd_data.condition{sub} == 1 | Todd_data.condition{sub} == 2 | Todd_data.condition{sub} == 5 | Todd_data.condition{sub} == 6);
    sh_trials = find(Todd_data.condition{sub} == 3 | Todd_data.condition{sub} == 4 | Todd_data.condition{sub} == 7 | Todd_data.condition{sub} == 8);
    strials(sl_trials) = -1;
    strials(sh_trials) = 1;
   Todd_data.social{sub} = strials';


    % assign 1 to conflict trials and -1 to noconflict trials  
    Todd_data.conflict{sub} = -1 * ones(96,1);
    Todd_data.conflict{sub}(find(Todd_data.condition{sub} > 2 & Todd_data.condition{sub} < 7)) = 1;
   

    % exact social stim values
  
    clear socmean socstd
    
    [dummy2 dummy] = xlsread(f, [f(1:4), '.txt'], 'BZ3:BZ103'); clear dummy2
    socstims = dummy([1:16 18:33 35:50 52:67 69:84 86:101]);
    

    %% Plot by temperature

    figure('Color', [1 1 1], 'Name', (subjects{sub}));
    
    subplot(3,3,1)
    scatter(Todd_data.temp{sub}, Todd_data.pain{sub});
    xlim([47.5 50.5])
        
    subplot(3,3,4)
    plot(find(Todd_data.temp{sub} == 48), Todd_data.pain{sub}(Todd_data.temp{sub} == 48), 'Color', [.4 .1 .1], 'LineWidth', 1.5); hold on;
    plot(find(Todd_data.temp{sub} == 49), Todd_data.pain{sub}(Todd_data.temp{sub} == 49), 'Color', [.6 .3 .3], 'LineWidth', 1.5); hold on;
    plot(find(Todd_data.temp{sub} == 50), Todd_data.pain{sub}(Todd_data.temp{sub} == 50), 'Color', [1 .4 .4],  'LineWidth', 1.5); hold on;
    legend('48C', '49C', '50C');

    subplot(3,3,7)
    bar([48,49,50], [nanmean(Todd_data.pain{sub}(Todd_data.temp{sub} == 48)), nanmean(Todd_data.pain{sub}(Todd_data.temp{sub} == 49)), nanmean(Todd_data.pain{sub}(Todd_data.temp{sub} == 50))]);
    xlim([47.5 50.5])

    %% Plot by condition
    
    clsltrials = find(Todd_data.condition{sub} == 1 | Todd_data.condition{sub} == 2);
    clshtrials = find(Todd_data.condition{sub} == 3 | Todd_data.condition{sub} == 4);
    chsltrials = find(Todd_data.condition{sub} == 5 | Todd_data.condition{sub} == 6);
    chshtrials = find(Todd_data.condition{sub} == 7 | Todd_data.condition{sub} == 8);

    %% EXPT Timecourse

%     figure('Color',[1 1 1],  'Name', (subjects{sub}));
    
    subplot(2,3,2)
    plot(clsltrials, Todd_data.exp{sub}(clsltrials), 'Color', cms(1,:), 'LineWidth',1.5); hold on; 
    plot(clshtrials, Todd_data.exp{sub}(clshtrials), '--', 'Color', cms(2,:), 'LineWidth',1.5); hold on; 
    plot(chsltrials, Todd_data.exp{sub}(chsltrials), '--', 'Color', cms(3,:), 'LineWidth',1.5); hold on; 
    plot(chshtrials, Todd_data.exp{sub}(chshtrials), 'Color', cms(4,:), 'LineWidth',1.5); hold on; 
    legend('CLSL', 'CLSH', 'CHSL', 'CHSH');
    title('Expectation ratings ALL trials','FontSize',16, 'FontName','Arial');
    
    % Barplots Expectations
    
    subplot(2,3,5)
    barplot_colored([Todd_data.exp{sub}(clsltrials), Todd_data.exp{sub}(clshtrials), Todd_data.exp{sub}(chsltrials), Todd_data.exp{sub}(chshtrials)])
    set(gca,'XTickLabel',{'CLSL', 'CLSH', 'CHSL', 'CHSH'});
    title('Expectation Ratings ALL trials','FontSize',16, 'FontName','Arial');
    colormap(cms)

    
    %% timecourse medium temp PAIN ratings by condition

%     figure('Color',[1 1 1],  'Name', (subjects{sub}));
    
    subplot(3,3,3)
    plot(clsltrials, Todd_data.pain{sub}(clsltrials), 'Color', cms(1,:), 'LineWidth',1.5); hold on; 
    plot(clshtrials, Todd_data.pain{sub}(clshtrials), '--', 'Color', cms(2,:), 'LineWidth',1.5); hold on; 
    plot(chsltrials, Todd_data.pain{sub}(chsltrials), '--', 'Color', cms(3,:), 'LineWidth',1.5); hold on; 
    plot(chshtrials, Todd_data.pain{sub}(chshtrials), 'Color', cms(4,:), 'LineWidth',1.5); hold on; 
    legend('CLSL', 'CLSH', 'CHSL', 'CHSH');
    title('PAIN ratings all heat trials','FontSize',16, 'FontName','Arial');

    subplot(3,3,6)
    barplot_colored([Todd_data.pain{sub}(Todd_data.condition{sub} == 2), Todd_data.pain{sub}(Todd_data.condition{sub} == 4), Todd_data.pain{sub}(Todd_data.condition{sub} == 5), Todd_data.pain{sub}(Todd_data.condition{sub} == 7)])
    set(gca,'XTickLabel',{'CLSL49', 'CLSH49', 'CHSL49', 'CHSH49'});
    title('PAIN Ratings medium heat','FontSize',16, 'FontName','Arial');
    colormap(cms)

    subplot(3,3,9)
    barplot_colored([Todd_data.pain{sub}(Todd_data.condition{sub} == 1), Todd_data.pain{sub}(Todd_data.condition{sub} == 3), Todd_data.pain{sub}(Todd_data.condition{sub} == 2), Todd_data.pain{sub}(Todd_data.condition{sub} == 4), Todd_data.pain{sub}(Todd_data.condition{sub} == 5), Todd_data.pain{sub}(Todd_data.condition{sub} == 7), Todd_data.pain{sub}(Todd_data.condition{sub} == 6), Todd_data.pain{sub}(Todd_data.condition{sub} == 8)])
    set(gca,'XTickLabel',{'CLSL48', 'CLSL49', 'CLSH48', 'CLSH49', 'CHSL49', 'CHSL50', 'CHSH49', 'CHSH50'});
    title('PAIN Ratings all','FontSize',16, 'FontName','Arial');
    colormap(cm)

    %%
    painmeans_todd = nanmean([Todd_data.pain{sub}(Todd_data.condition{sub} == 1), Todd_data.pain{sub}(Todd_data.condition{sub} == 3), Todd_data.pain{sub}(Todd_data.condition{sub} == 2), Todd_data.pain{sub}(Todd_data.condition{sub} == 4), Todd_data.pain{sub}(Todd_data.condition{sub} == 5), Todd_data.pain{sub}(Todd_data.condition{sub} == 7), Todd_data.pain{sub}(Todd_data.condition{sub} == 6), Todd_data.pain{sub}(Todd_data.condition{sub} == 8)]);
    expemeans_todd = nanmean([Todd_data.exp{sub}(clsltrials), Todd_data.exp{sub}(clshtrials), Todd_data.exp{sub}(chsltrials), Todd_data.exp{sub}(chshtrials)]);
    pain_todd_temp = [nanmean(Todd_data.pain{sub}(Todd_data.temp{sub} == 48)), nanmean(Todd_data.pain{sub}(Todd_data.temp{sub} == 49)), nanmean(Todd_data.pain{sub}(Todd_data.temp{sub} == 50))];

    
%%
input('Press any key to continue to grand average')


%% load other people's data
cd ('/Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Behavioral/Results')
load SCEBLmri_data.mat


%% Temp effects
figure('Color', [1 1 1])
barplot_colored(pain_all_temp, 'within'); hold on
plot([1 2 3], pain_all_temp, 'ko'); hold on
set(gca,'XTickLabel',{'48C', '49C', '50C'});
title('Pain ratings by temperatures','FontSize',16, 'FontName','Arial');
colormap(cmtemp)

% add here Todd's pain ratings in red
plot([1 2 3], pain_todd_temp, 'r*'); hold on

%% Overall pain sensitiviy (average pain ratings)

figure('Color', [1 1 1])

plot(mean(painmeans_all,2), ones(30,1), 'ok'); hold on
plot(mean(painmeans_todd), 1, '*r')
xlim([0 100])


%% Expectations and medium pain

figure('Color', [1 1 1]) 
subplot(2,1,1)
barplot_colored(expemeans_all, 'within'); hold on
plot([1:4], expemeans_todd, 'ow')
set(gca,'XTickLabel',{'CLSL', 'CLSH', 'CHSL', 'CHSH'});
title('Expectation Ratings','FontSize',16, 'FontName','Arial');
colormap(cms)

subplot(2,1,2)
barplot_colored(painmeans_all(:,3:6), 'within'); hold on
plot([1:4], painmeans_todd(3:6), 'ow')
set(gca,'XTickLabel',{'T49 CLSL', 'T49 CLSH', 'T49 CHSL', 'T49 CHSH'});
title('Pain Ratings medium temperatures','FontSize',16, 'FontName','Arial');
colormap(cms)

%% All pain

figure('Color', [1 1 1]) 
barplot_colored(painmeans_all, 'within'); hold on
plot([1:8], painmeans_todd, 'ow')
set(gca,'XTickLabel',{'T48 CLSL', 'T48 CLSH', 'T49 CLSL', 'T49 CLSH', 'T49 CHSL', 'T49 CHSH', 'T50 CHSL', 'T50 CHSH'});
title('Pain Ratings ALL trials and temperatures','FontSize',16, 'FontName','Arial');
colormap(cm)


%% social influence effect for all temperatures

figure('Color', [1 1 1]) 

subplot(1,2,1)
barplot_colored([painmeans_all(:,1), mean([painmeans_all(:,3), painmeans_all(:,5)], 2), painmeans_all(:,7)], 'within'); hold on
plot(1:3, [painmeans_todd(1), mean([painmeans_todd(3), painmeans_todd(5)], 2), painmeans_todd(:,7)], 'ow')
ylim([0 50])

subplot(1,2,2)
barplot_colored([painmeans_all(:,2), mean([painmeans_all(:,4), painmeans_all(:,6)], 2), painmeans_all(:,8)], 'within'); hold on
plot(1:3, [painmeans_todd(2), mean([painmeans_todd(4), painmeans_todd(6)], 2), painmeans_todd(:,8)], 'ow')
colormap(cm)
ylim([0 50])


%% Mindfulness pre vs post

figure('Color', [1 1 1]) ;
bar(1:2, [nanmean(Todd_data.pain{1}(1:48)), nanmean(Todd_data.pain{1}(49:96))])





%% Time course learning on average

% One figure for time course of expectation and pain ratings

figure('Color', [1 1 1]) ;
subplot(2,1,1)
plot(1:24, exp_tc.clsl, 'Color', cms(1,:), 'LineWidth', 1.5); hold on
plot(1:24, exp_tc.clsh, 'Color', cms(2,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
plot(1:24, exp_tc.chsl, 'Color', cms(3,:), 'LineWidth', 1.5); hold on
plot(1:24, exp_tc.chsh, 'Color', cms(4,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
title('Expectation ratings by time', 'FontSize', 16, 'FontName', 'Arial');
xlim([0 25])
legend('clsl', 'clsh', 'chsl', 'chsh');

subplot(2,1,2)
plot(1:12, pain49_tc.clsl, 'Color', cms(1,:), 'LineWidth', 1.5); hold on
plot(1:12, pain49_tc.clsh, 'Color', cms(2,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
plot(1:12, pain49_tc.chsl, 'Color', cms(3,:), 'LineWidth', 1.5); hold on
plot(1:12, pain49_tc.chsh, 'Color', cms(4,:), 'LineWidth', 1.5, 'LineStyle', '--'); hold on
title('Pain ratings medium temp by time', 'FontSize', 16, 'FontName', 'Arial');
xlim([0.5 12.5])
legend('clsl', 'clsh', 'chsl', 'chsh');



%% regression lines for relationship between expectation and pain ratings

figure('Color', [1 1 1])
h = line_plot_multisubject(Todd_data.exp, Todd_data.pain); title('Multi-subject scatterplot: expectation ratings - pain ratings','FontSize',16, 'FontName', 'Arial', 'FontWeight', 'Bold')
xlabel('Expectation', 'FontSize', 14, 'FontName','Arial');
ylabel('Pain', 'FontSize', 14, 'FontName','Arial');


cd /Users/Leonie/Documents/BOULDER/PROJECTS/6_SCEBL1_fMRI/Behavioral/Results
save SCEBLmri_data_Todd Todd_data


%% load other people's data

load SCEBLmri_data.mat

%% Temp effects
figure('Color', [1 1 1])
barplot_colored(pain_all_temp, 'within'); hold on
plot([1 2 3], pain_all_temp, 'ko'); hold on
set(gca,'XTickLabel',{'48C', '49C', '50C'});
title('Pain ratings by temperatures','FontSize',16, 'FontName','Arial');
colormap(cmtemp)

%% add here Todd's pain ratings in red

plot([1 2 3], [20 40 60], 'r*'); hold on

% output mean and 95% CI values for temperatures
temp_means = mean(pain_all_temp)
temp_ste = ste(pain_all_temp)
temp_95ci_l = temp_means - 2.001 * temp_ste
temp_95ci_h = temp_means + 2.001 * temp_ste






