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

[p2stats] = processTABLE(p2 , 2);

[p3stats] = processTABLE(p3 , 3); 

p1stats = array2table(p1stats); 
p1stats.Properties.RowNames = {'P value', 'Bootstrap KS', 'Bootstrap W', 'H'}; 
p1stats.Properties.VariableNames = {'L12', 'L13', 'L23', 'R12', 'R13', 'R23'};  

p2stats = array2table(p2stats); 
p2stats.Properties.RowNames = {'P value', 'Bootstrap KS', 'Bootstrap W', 'H'}; 
p2stats.Properties.VariableNames = {'L12', 'L13', 'L23', 'R12', 'R13', 'R23'};  

p3stats = array2table(p3stats); 
p3stats.Properties.RowNames = {'P value', 'Bootstrap KS', 'Bootstrap W', 'H'}; 
p3stats.Properties.VariableNames = {'L12', 'L13', 'L23', 'R12', 'R13', 'R23'};  


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

t = tiledlayout(2, 3); 
% Plot Kernal Density using probability distribution function (default for
% function)
[f_L1,xi_L1,~] = ksdensity(left_1);
[f_L2,xi_L2,~] = ksdensity(left_2);
[f_L3,xi_L3,~] = ksdensity(left_3);

[f_R1,xi_R1,~] = ksdensity(right_1);
[f_R2,xi_R2,~] = ksdensity(right_2);
[f_R3,xi_R3,~] = ksdensity(right_3);

% LEFT and RIGHT - Run 1 vs 2
nexttile(1)
plot(xi_L1,f_L1);
hold on
plot(xi_L2,f_L2);
title(['p value ', num2str(LEFT_r12)])
ylabel('Probability density')
xlabel('Log normal power')
legend([pLabel, ' Left R1'],[pLabel, ' Left R2'])

nexttile(4);
plot(xi_R1,f_R1);
hold on
plot(xi_R2,f_R2);
title(['p value ', num2str(RIGHT_r12)])
ylabel('Probability density')
xlabel('Log normal power')
legend([pLabel, ' Right R1'],[pLabel, ' Right R2'])



% LEFT and RIGHT - Run 1 vs 3
nexttile(2)
plot(xi_L1,f_L1);
hold on
plot(xi_L3,f_L3);
title(['p value ', num2str(LEFT_r13)])
ylabel('Probability density')
xlabel('Log normal power')
legend([pLabel, ' Left R1'],[pLabel, ' Left R3'])

nexttile(5)
plot(xi_R1,f_R1);
hold on
plot(xi_R3,f_R3);
title(['p value ', num2str(RIGHT_r13)])
ylabel('Probability density')
xlabel('Log normal power')
legend([pLabel, ' Right R1'],[pLabel, ' Right R3'])


% LEFT and RIGHT - Run 2 vs 3
nexttile(3)
plot(xi_L2,f_L2);
hold on
plot(xi_L3,f_L3);
title(['p value ', num2str(LEFT_r23)])
ylabel('Probability density')
xlabel('Log normal power')
legend([pLabel, ' Left R2'],[pLabel, ' Left R3'])

nexttile(6)
plot(xi_R2,f_R2);
hold on
plot(xi_R3,f_R3);
title(['p value ', num2str(RIGHT_r23)])
ylabel('Probability density')
xlabel('Log normal power')
legend([pLabel, ' Right R2'],[pLabel, ' Right R3'])

outStats = [LEFT_r12 , LEFT_r13 , LEFT_r23 , RIGHT_r12 , RIGHT_r13 , RIGHT_r23];

t.TileSpacing = 'compact';
t.Padding = 'compact';




figure
t2 = tiledlayout(2, 3);

%LEFT 
[~, ~, KS12] = kstest2(left_R1, left_R2);
[~, ~, KS13] = kstest2(left_R1, left_R3);
[~, ~, KS23] = kstest2(left_R2, left_R3);

[~, ~, WS12] = ranksum(left_R1, left_R2);
[~, ~, WS13] = ranksum(left_R1, left_R3);
[~, ~, WS23] = ranksum(left_R2, left_R3);

r1 = left_R1;
r2 = left_R2;
r3 = left_R3;

%RUN 1 VS 2 LEFT 
ksvals = zeros(1000,1);
wsvals = zeros(1000,1);
for i = 1:1000
    x = round(length(r1)/10); %10% of data
    rindex = randsample(length(r1), x); %get random once and use it for both
    a = r1(rindex);
    b = r2(rindex);
    [~,~,wsstats] = ranksum(a,b);
    [~,~,ksstats] = kstest2(a,b);
    ksvals(i) = ksstats;
    wsvals(i) = wsstats.zval;
end
P12l = sum(ksvals < KS12) / length(ksvals);
W12l = sum(wsvals < WS12.zval) / length(wsvals);

if P12l && W12l > 0.05
    H12l = 0;
else if P12l && W12l < 0.05
        H12l = 1;
    else
        H12l = 7;
    end
end 

nexttile(1)
histogram(ksvals)
hold on 
xline(KS12)
title(["L12", P12l, pLabel])


%RUN 1 VS 3 LEFT 
ksvals = zeros(1000,1);
wsvals = zeros(1000,1);
for i = 1:1000
    x = round(length(r1)/10); %10% of data
    rindex = randsample(length(r1), x); %get random once and use it for both
    a = r1(rindex);
    b = r3(rindex);
    [~,~,wsstats] = ranksum(a,b);
    [~,~,ksstats] = kstest2(a,b);
    ksvals(i) = ksstats;
    wsvals(i) = wsstats.zval;
end
P13l = sum(ksvals < KS13) / length(ksvals);
W13l = sum(wsvals < WS13.zval) / length(wsvals);

if P13l && W13l > 0.05
    H13l = 0;
else if P13l && W13l < 0.05
        H13l = 1;
    else
        H13l = 7;
    end
end 

nexttile(2)
histogram(ksvals)
hold on 
xline(KS13)
title(["L13", P13l, pLabel])

%RUN 2 VS 3 LEFT 
ksvals = zeros(1000,1);
wsvals = zeros(1000,1);
for i = 1:1000
    x = round(length(r2)/10); %10% of data
    rindex = randsample(length(r2), x); %get random once and use it for both
    a = r2(rindex);
    b = r3(rindex);
    [~,~,wsstats] = ranksum(a,b);
    [~,~,ksstats] = kstest2(a,b);
    ksvals(i) = ksstats;
    wsvals(i) = wsstats.zval;
end
P23l = sum(ksvals < KS23) / length(ksvals);
W23l = sum(wsvals < WS23.zval) / length(wsvals);

if P23l > 0.05 && W23l > 0.05
    H23l = 0;
else if P23l < 0.05 && W23l < 0.05
        H23l = 1;
    else
        H23l = 7;
    end
end 

nexttile(3)
histogram(ksvals)
hold on 
xline(KS23)
title(["L23", P23l, pLabel])

%RIGHT SIDE  
[~, ~, KS12] = kstest2(right_R1, right_R2);
[~, ~, KS13] = kstest2(right_R1, right_R3);
[~, ~, KS23] = kstest2(right_R2, right_R3);

[~, ~, WS12] = ranksum(right_R1, right_R2);
[~, ~, WS13] = ranksum(right_R1, right_R3);
[~, ~, WS23] = ranksum(right_R2, right_R3);

r1 = right_R1;
r2 = right_R2;
r3 = right_R3;

%RUN 1 VS 2 LEFT 
ksvals = zeros(1000,1);
wsvals = zeros(1000,1);
for i = 1:1000
    x = round(length(r1)/10); %10% of data
    rindex = randsample(length(r1), x); %get random once and use it for both
    a = r1(rindex);
    b = r2(rindex);
    [~,~,wsstats] = ranksum(a,b);
    [~,~,ksstats] = kstest2(a,b);
    ksvals(i) = ksstats;
    wsvals(i) = wsstats.zval;
end
P12r = sum(ksvals < KS12) / length(ksvals);
W12r = sum(wsvals < WS12.zval) / length(wsvals);

if P12r && W12r > 0.05
    H12r = 0;
else if P12r && W12r < 0.05
        H12r = 1;
    else
        H12r = 7;
    end
end

nexttile(4)
histogram(ksvals)
hold on 
xline(KS12)
title(["R12", P12r, pLabel])

%RUN 1 VS 3 LEFT 
ksvals = zeros(1000,1);
wsvals = zeros(1000,1);
for i = 1:1000
    x = round(length(r1)/10); %10% of data
    rindex = randsample(length(r1), x); %get random once and use it for both
    a = r1(rindex);
    b = r3(rindex);
    [~,~,wsstats] = ranksum(a,b);
    [~,~,ksstats] = kstest2(a,b);
    ksvals(i) = ksstats;
    wsvals(i) = wsstats.zval;
end
P13r = sum(ksvals < KS13) / length(ksvals);
W13r = sum(wsvals < WS13.zval) / length(wsvals);

if P13r && W13r > 0.05
    H13r = 0;
else if P13r && W13r < 0.05
        H13r = 1;
    else
        H13r = 7;
    end
end 

nexttile(5)
histogram(ksvals)
hold on 
xline(KS13)
title(["L13", P13r, pLabel])

%RUN 2 VS 3 LEFT 
ksvals = zeros(1000,1);
wsvals = zeros(1000,1);
for i = 1:1000
    x = round(length(r2)/10); %10% of data
    rindex = randsample(length(r2), x); %get random once and use it for both
    a = r2(rindex);
    b = r3(rindex);
    [~,~,wsstats] = ranksum(a,b);
    [~,~,ksstats] = kstest2(a,b);
    ksvals(i) = ksstats;
    wsvals(i) = wsstats.zval;
end
P23r = sum(ksvals < KS23) / length(ksvals);
W23r = sum(wsvals < WS23.zval) / length(wsvals);

if P23r > 0.05 && W23r > 0.05
    H23r = 0;
else if P23r < 0.05 && W23r < 0.05
        H23r = 1;
    else
        H23r = 7;
    end
end 

nexttile(6)
histogram(ksvals)
hold on 
xline(KS23)
title(["L23", P23r, pLabel])

outKS = [LEFT_r12 , LEFT_r13 , LEFT_r23 , RIGHT_r12 , RIGHT_r13 , RIGHT_r23];
outP = [P12l, P13l, P23l, P12r, P13r, P23r]; 
outW = [W12l, W13l, W23l, W12r, W13r, W23r];  
outH = [H12l, H13l, H23l, H12r, H13r, H23r]; 
outStats = [outKS; outP; outW; outH]; 
end 