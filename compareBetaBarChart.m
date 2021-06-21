function [maxBetas] = compareBetaBarChart(startindex, leng, js, channel)
%start index may have to be manually adjusted if junk trials happen
%js is the json file to be read
%channel is the cell array to get the data for

ifcounter = 0;          %used to put each trial in a new column
y = zeros(4096,3);      %length fits output from pspectruum, 3 columns for each trial
locallength = length(channel);
maxBetas = zeros(6, 1);
ystd = zeros(6, 1);
yvar = zeros(6, 1);
for m = 1:locallength
    for i = startindex:leng
        if strcmp(js.LfpMontageTimeDomain(i).Channel,channel{m}) 
            ifcounter = ifcounter + 1;
            t = js.LfpMontageTimeDomain(i).TimeDomainData;
            [p,f] = pspectrum(t, 250, 'FrequencyLimits', [0 100]); %250 comes from json file itself
            y(:,ifcounter) = p;
        end
    end 
    ym = mean(y, 2); %average across runs (by rows)
    fbeta = f > 12 & f < 33; %logical to only find beta region 
    maxBetas(m) = max(ym(fbeta)); %find max within ym within beta 
    ystd(m) = std(ym(fbeta));
    yvar(m) = var(ym(fbeta));
    ifcounter = 0; 
end
%label the max

figure
bar(maxBetas)
set(gca,'xticklabel',{'03', '13', '02', '01', '12', '23'})
hold on
errorbar(1:6, maxBetas, ystd)
title("Comparing max Betas (std)")

figure 
bar(maxBetas)
set(gca,'xticklabel',{'03', '13', '02', '01', '12', '23'})
hold on
errorbar(1:6, maxBetas, yvar)
title("Comparing max Betas (var)")
end