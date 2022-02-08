%% Change directory to the folder that contains the Patient CSV files
%JAT machine
%cd('I:\01_Coding_Datasets\PERCEPT MDT colab\BrainSense Survey Stats')
%SL machine
cd('C:\Users\sydne\Documents\github\perceive\patientData\abstract_data')

%% Read in the CSV files for each patient
p1 = readtable('brainSenseSurvey_P1.csv');
p2 = readtable('brainSenseSurvey_P2.csv');
p3 = readtable('brainSenseSurvey_P3.csv');

%% Compute the stats for each patient and create plots
close all

[p1stats] = processTABLE(p1, 1);

%[p2stats] = processTABLE(p2 , 2);

%[p3stats] = processTABLE(p3 , 3);

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

t = tiledlayout(3,2);
% Plot Kernal Density using probability distribution function (default for
% function)
[f_L1,xi_L1,~] = ksdensity(left_1);
[f_L2,xi_L2,~] = ksdensity(left_2);
[f_L3,xi_L3,~] = ksdensity(left_3);

[f_R1,xi_R1,~] = ksdensity(right_1);
[f_R2,xi_R2,~] = ksdensity(right_2);
[f_R3,xi_R3,~] = ksdensity(right_3);

% LEFT and RIGHT - Run 1 vs 2
nexttile
plot(xi_L1,f_L1);
hold on
plot(xi_L2,f_L2);
%title(['p value ', num2str(LEFT_r12, 3)])
title(['left Run 1 vs Run 2, p = ', num2str(round(LEFT_r12, 3))])
ylabel('Probability density')
xlabel('Log normal power')
%legend([pLabel, ' Left R1'],[pLabel, ' Left R2'])
legend([pLabel, 'Run 1'],[pLabel, ' Run 2'])

nexttile
plot(xi_R1,f_R1);
hold on
plot(xi_R2,f_R2);
%title(['p value ', num2str(RIGHT_r12, 3)])
title(['right Run 1 vs Run 2, p = ', num2str(round(RIGHT_r12, 3))])
ylabel('Probability density')
xlabel('Log normal power')
%legend([pLabel, ' Left R1'],[pLabel, ' Left R2'])
legend([pLabel, ' Run 1'],[pLabel, ' Run 2'])


% LEFT and RIGHT - Run 1 vs 3
nexttile
plot(xi_L1,f_L1);
hold on
plot(xi_L3,f_L3);
%title(['p value ', num2str(LEFT_r13, 3)])
title(['left Run 1 vs Run 3, p = ', num2str(round(LEFT_r13, 3))])
ylabel('Probability density')
xlabel('Log normal power')
%legend([pLabel, ' Left R1'],[pLabel, ' Left R2'])
legend([pLabel, ' Run 1'],[pLabel, ' Run 3'])


nexttile
plot(xi_R1,f_R1);
hold on
plot(xi_R3,f_R3);
%title(['p value ', num2str(RIGHT_r13, 3)])
title(['right Run 1 vs Run 3, p = ', num2str(round(RIGHT_r13, 3))])
ylabel('Probability density')
xlabel('Log normal power')
%legend([pLabel, ' Left R1'],[pLabel, ' Left R2'])
legend([pLabel, ' Run 1'],[pLabel, ' Run 3'])


% LEFT and RIGHT - Run 2 vs 3
nexttile
plot(xi_L2,f_L2);
hold on
plot(xi_L3,f_L3);
%title(['p value ', num2str(LEFT_r23, 3)])
title(['left Run 2 vs Run 3, p = ', num2str(round(LEFT_r23, 3))])
ylabel('Probability density')
xlabel('Log normal power')
%legend([pLabel, ' Left R1'],[pLabel, ' Left R2'])
legend([pLabel, ' Run 2'],[pLabel, ' Run 3'])


nexttile
plot(xi_R2,f_R2);
hold on
plot(xi_R3,f_R3);
%title(['p value ', num2str(RIGHT_r23, 3)])
title(['right Run 2 vs Run 3, p = ', num2str(round(RIGHT_r23, 3))])
ylabel('Probability density')
xlabel('Log normal power')
%legend([pLabel, ' Left R1'],[pLabel, ' Left R2'])
legend([pLabel, ' Run 2'],[pLabel, ' Run 3'])


outStats = [LEFT_r12 , LEFT_r13 , LEFT_r23 , RIGHT_r12 , RIGHT_r13 , RIGHT_r23];

t.TileSpacing = 'compact';
t.Padding = 'compact'; 

% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.25, 0.98, 'Patient 1 Left Hemisphere') 

sleft = sgtitle("Patient 1 Left Hemisphere"); 
sleft.HorizontalAlignment = "right"; 
hold on
sright = sgtitle("Patient 1 Right Hemisphere"); 
sright.HorizontalAlignment = "left"; 


end