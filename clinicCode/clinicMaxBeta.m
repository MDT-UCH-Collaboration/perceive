%NEXT STEPS
%create some sort of function (probably using the diff duration array) 
%to figure out if and where junk trials happened 

%clean out the random code towards the end, only keeping it bc it might 
%be useful to the point above 

%create a smarter system for the compareBetaBarChart and allFreqBarChart
%based on it its double or single battery 
%i.e don't hard code "channelsRight" or "channelsLeft" 

%impt note: the channels are in a STABLE arrangement! in some places 
%(althought I think only compBetaBarChart, the order of channels 
%is hard coded, being 03, 13, 02, etc 
%if the stable thing gets changed then be mindful about where the order 
%is hard coded
%% loads in json file
%laods in file 
cd('C:\Users\sydne\Documents\MATLAB\ThompsonLab\Patient1_0524')
jsonFiles = 'Report_Json_Session_Report_20210524T101039.json';
js = jsondecode(fileread(jsonFiles));
%sets it back to the path where all the other functions are 
cd('C:\Users\sydne\Documents\MATLAB\ThompsonLab')

%declares "global" variables AND determines if .json comes from 
%single or double battery B)
channels = unique({js.LfpMontageTimeDomain.Channel}, 'stable');
leng = numel({js.LfpMontageTimeDomain.Channel});
sides = {'LEFT', 'RIGHT'};
if any((contains(channels, sides{1}))) && any((contains(channels, sides{2}))) 
        doubleBattery = false; 
        %json has both left and right data 
        channelsLeft = channels(contains(channels, sides{1}));
        channelsRight = channels(contains(channels, sides{2}));
    else 
        doubleBattery = true; 
end 

%% create table showing specfics for each run and channel 

maxMeanBetaAllRuns(1, js, leng, channelsRight)
%% run to get bar chart comparing betas w std and var 
testvar = compareBetaBarChart(7, leng, js, channels);
%% run to get full graph
%{} gets the thing from the actual cell array
%() gets cell ARRAY!

figure;

if doubleBattery == false
    %this should be true if the json returns both left and right data,
    %single battery
    for i = 1:2
        subplot(2, 1, i)
        channels2 = channels(contains(channels,sides{i}));
        colors = 'krbgmc';
        maxValues = zeros(length(channels2), 1);
        for c=1:length(channels2)
            maxValues(c) = tempAvgPlot(1, leng, js, channels2{c}, colors(c));
            hold on;
        end
        
        highestBetas = {' ', ' '};
        [M, I] = max(maxValues);
        highestBeta = channels2{I};
        highestBetas{i} = highestBeta;
        
        %add all the plot info
        xline(13)
        xline(30)
        
        legendChan = cell(6);
        oldchar = {'LEFT', 'RIGHT', '0', '1', '2', '3', '4', '5', '_', 'AND'}; 
        newchar = {''}; 
        oldnum = {'ZERO', 'ONE', 'TWO', 'THREE'}; 
        newnum = {'0', '1', '2', '3'};
        for b = 1:length(channels2)
            legendChan{b} = replace(channels2{b}, oldchar, newchar);
            legendChan{b} = replace(legendChan{b}, oldnum, newnum);
        end
        legend(legendChan{1}, legendChan{2}, legendChan{3}, legendChan{4}, legendChan{5}, legendChan{6})
        
        title(["Freq vs. Power (db)", sides{i}, highestBeta])
        xlim([0 60])
        xlabel("Frequency")
        ylabel("Power")
        
        disp(['The highest beta comes from ', highestBeta])
    end
    
    
else
    %should be true if a double battery
    x = 0;
    %leng = numel({js.LfpMontageTimeDomain.Channel});
    %channels2 = channels(contains(channels,sides{i}));
    
    
    colors = 'krbgmc';
    maxValues = zeros(length(channels), 1);
    for c=1:length(channels)
        maxValues(c) = tempAvgPlot(7, leng, js, channels{c}, colors(c));
        hold on;
    end 
    
    [M, I] = max(maxValues);
    highestBeta = channels{I};
    
    %add all the plot info
    
    xline(13);
    xline(30);
    
    %these are the order presented in channels
    
    legendChan = {};
    for i = 1:length(channels)
        legendChan{i} = replace(channels{i}, '_', ' ');
    end
    legend(legendChan{1}, legendChan{2}, legendChan{3}, legendChan{4}, legendChan{5}, legendChan{6}, "Beta lb", "Beta ub")
    
    %legend("03", "13", "02", "01", "12", "23")
    titlehb = replace(highestBeta, '_', ' ');
    title(["Freq vs. Power (db)", legendChan{I}])
    xlabel("Frequency")
    ylabel("Power")
    
    disp(['The highest beta comes from ', highestBeta])
end
%% Run to Create bar chart summary for all frequency bands of max beta channel
%would have to go into if statement
allFreqBarChart(7, leng, js, highestBeta) 
%% Direct compare two contact pairs
close;
channelA = channels{2};
channelB = channels{5};
directCompare(7, leng, js, channelA, channelB)

xline(13)
xline(30)
%these are the order presented in channels
legend(channelA, channelB, "Beta lb", "Beta ub")
title("Freq vs. Power (db)")
xlabel("Frequency")
ylabel("Power")
%% using function to compare specific trials
%close;
%numbers should be six apart spanning the trials
figure;
nonAvgTrialCompare(leng, js, 1, 13, 13+12)
xlabel("Frequency (Hz)")
ylabel("Power Spectrum (dB)")
legend("1", "13", "25")
figure;
nonAvgTrialCompare(leng, js, 8, 20, 32)
xlabel("Frequency (Hz)")
ylabel("Power Spectrum (dB)")
legend("8", "20", "32")
%% Not sure if this is relevant any more? 
%Generates differences between trials
testtime = {js.LfpMontageTimeDomain.FirstPacketDateTime};
cleantime = {}; %initialize new cell array for clean times
leng = numel(testtime); %have to use numel here bc length sucks?
for i = 1:leng
    bettertime = replace(testtime{i}, {'T', 'Z'},  {' ', ''});
    cleantime{i} = datetime(bettertime) - hours(6); %converts string to datetime
    %fix by moving the time back six hours to make it correct
end

%now find differentials between the times
diff = {}; %q is this still a cell array??
for i = 1:(leng-1)
    diff{i} = cleantime{i+1} - cleantime{i};
end
%the first entry is the difference between cells 2 and 1,
%second entry is diff between 3 and 2, etc
