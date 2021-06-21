function [maxBeta] = tempAvgPlot(startindex, leng, js, channelnum, color) 
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
f = f';  %transposing both to fit size requirements set by shadedErrorBar
y = y';

smy = smoothdata(y, 2, "sgolay", 250); %y is now smooth
%shadedErrorBar(f, smy, {@mean, @std}, 'lineprops', ['*', color], 'transparent', false);
ym = mean(y, 1); %takes the mean of y  
testym = ym'; 
plot(f, ym, 'color', color);
%label the max 
isBeta = f > 12 & f < 33; %creates logical matrix within beta region
check = (ym(isBeta)); 
maxBeta = max(ym(isBeta)); %finds power max within the beta region


end 