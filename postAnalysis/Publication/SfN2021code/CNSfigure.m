%% Change directory to the folder that contains the Patient CSV files
%cd('I:\01_Coding_Datasets\PERCEPT MDT colab\BrainSense Survey Stats')
cd('C:\Users\sydne\Documents\github\perceive\patientData\abstract_data')
%% Read in the CSV files for each patient
p1 = readtable('brainSenseSurvey_P1.csv');
p2 = readtable('brainSenseSurvey_P2.csv');
p3 = readtable('brainSenseSurvey_P3.csv');

%% Compute the stats for each patient and create plots
close all

[p1stats] = processTABLE(p1, 1);

% [p2stats] = processTABLE(p2 , 2);
% 
% [p3stats] = processTABLE(p3 , 3);

%%

function [newTab] = transformTAB(oldTab)

oldTab.side = string(oldTab.side);
% Convert 'Side' from 1 = left and 2 = right [hemisphere]
oldTab.side(ismember(oldTab.side,"1"),:) = "L";
oldTab.side(ismember(oldTab.side,"2"),:) = "R";
% Transform power to absolute log 10 power
oldTab.betapwr = abs(log10(oldTab.betapwr));


newTab = oldTab;

end




function [outStats] = processTABLE(inP , pNUM)

figure;

[newTab] = transformTAB(inP);

pLabel = ['P' , num2str(pNUM)];

% Extract hemisphere and trial run number 
left_R1 = newTab.betapwr(contains(newTab.side,"L") & ismember(newTab.run, 1),:);
left_R2 = newTab.betapwr(contains(newTab.side,"L") & ismember(newTab.run, 2),:);
left_R3 = newTab.betapwr(contains(newTab.side,"L") & ismember(newTab.run, 3),:);
right_R1 = newTab.betapwr(contains(newTab.side,"R") & ismember(newTab.run, 1),:);
right_R2 = newTab.betapwr(contains(newTab.side,"R") & ismember(newTab.run, 2),:);
right_R3 = newTab.betapwr(contains(newTab.side,"R") & ismember(newTab.run, 3),:);

% Downsample from 250 Hz to 25 Hz
left_1 = downsample(left_R1,10);
left_2 = downsample(left_R2,10);
left_3 = downsample(left_R3,10);
% Downsample from 250 Hz to 25 Hz
right_1 = downsample(right_R1,10);
right_2 = downsample(right_R2,10);
right_3 = downsample(right_R3,10);

% Run kolmogorov smirnov 2 distribution test
% example name - LEFT hemisphere run 1 vs run 2 [LEFT_r12]
[~,LEFT_r12,~] = kstest2(left_1,left_2);
[~,LEFT_r13,~] = kstest2(left_1,left_3);
[~,LEFT_r23,~] = kstest2(left_2,left_3);
% Run Kolgomorov-Smirnov 2 distribution test
[~,RIGHT_r12,~] = kstest2(right_1,right_2);
[~,RIGHT_r13,~] = kstest2(right_1,right_3);
[~,RIGHT_r23,~] = kstest2(right_2,right_3);

t = tiledlayout(3,1); 
% Plot Kernal Density using probability distribution function (default for
% function)
[f_L1,xi_L1,~] = ksdensity(left_1);
[f_L2,xi_L2,~] = ksdensity(left_2);
[f_L3,xi_L3,~] = ksdensity(left_3);

% [f_R1,xi_R1,~] = ksdensity(right_1);
% [f_R2,xi_R2,~] = ksdensity(right_2);
% [f_R3,xi_R3,~] = ksdensity(right_3);

% LEFT and RIGHT - Run 1 vs 2
nexttile(1)
plot(xi_L1,f_L1);
hold on
plot(xi_L2,f_L2);
title(['p value = ', num2str(LEFT_r12)])
ylabel('Probability density')
xlabel('Log normal power')
legend(['Run 1 '],['Run 2'])

% nexttile(4);
% plot(xi_R1,f_R1);
% hold on
% plot(xi_R2,f_R2);
% title(['p value ', num2str(RIGHT_r12)])
% ylabel('Probability density')
% xlabel('Log normal power')
% legend([pLabel, ' Right R1'],[pLabel, ' Right R2'])



% LEFT and RIGHT - Run 1 vs 3
nexttile(2)
plot(xi_L1,f_L1);
hold on
plot(xi_L3,f_L3);
title(['p value = ', num2str(LEFT_r13)])
ylabel('Probability density')
xlabel('Log normal power')
legend(['Run 1'],['Run 3'])

% nexttile(5)
% plot(xi_R1,f_R1);
% hold on
% plot(xi_R3,f_R3);
% title(['p value ', num2str(RIGHT_r13)])
% ylabel('Probability density')
% xlabel('Log normal power')
% legend([pLabel, ' Right R1'],[pLabel, ' Right R3'])


% LEFT and RIGHT - Run 2 vs 3
nexttile(3)
plot(xi_L2,f_L2);
hold on
plot(xi_L3,f_L3);
title(['p value = ', num2str(LEFT_r23)])
ylabel('Probability density')
xlabel('Log normal power')
legend(['Run 2'],['Run 3'])
% 
% nexttile(6)
% plot(xi_R2,f_R2);
% hold on
% plot(xi_R3,f_R3);
% title(['p value ', num2str(RIGHT_r23)])
% ylabel('Probability density')
% xlabel('Log normal power')
% legend([pLabel, ' Right R2'],[pLabel, ' Right R3']) 

sgt = sgtitle(["Run comparision for LFP bipolar assessment"; "Patient 1 Left Hemisphere"]); 

outStats = [LEFT_r12 , LEFT_r13 , LEFT_r23 , RIGHT_r12 , RIGHT_r13 , RIGHT_r23];

t.TileSpacing = 'compact';
t.Padding = 'compact'; 
end