function [] = allFreqBarChart(startindex, leng, js, channelnum)
%start index may have to be manually adjusted if junk trials happen
%js is the json file to be read
%channelnum is the specific channel in the lfpmontagetimedomain to be read
%color creates the color of the line

ifcounter = 0;          %used to put each trial in a new column
y = zeros(4096,3);      %length fits output from pspectruum, 3 columns for each trial
for i = startindex:leng
    if strcmp(js.LfpMontageTimeDomain(i).Channel,channelnum) %might have to address differently
        ifcounter = ifcounter + 1;
        t = js.LfpMontageTimeDomain(i).TimeDomainData;
        [p,f] = pspectrum(t, 250, 'FrequencyLimits', [0 100]); %250 comes from json file itself
        y(:,ifcounter) = p;
    end
end

figure;
ym = mean(y, 2); %takes the mean of y by rows
% isBeta = f > 12 & f < 33; %creates logical matrix within beta region
% maxBeta = max(ym(isBeta)); %finds power max within the beta region
edges = [0 4; 5 8; 9 12; 13 30; 31 100];  %how many numbers between 0 to 4 exit in ym
%just use a loop, index f by the thresholds for bands
ymean = zeros(5, 1);
ystd = zeros(5, 1);
yvar = zeros(5, 1);
%ranges = categorical({'delta', 'theta', 'alpha', 'beta', 'gamma'});
for b = 1:size(edges, 1)
    strf = edges(b, 1);
    stpf = edges(b, 2);
    fnew = f > strf & f < stpf;
    ymean(b) = mean(ym(fnew));
    ystd(b) = std(ym(fnew));
    yvar(b) = var(ym(fnew));
end
figure
bar(ymean)
set(gca,'xticklabel',{'delta','theta','alpha','beta','gamma'})
hold on
errorbar(1:5, ymean, ystd)
titlechannel = string(replace(channelnum, '_', ' '));
title(["All frequency info for", titlechannel])
end 