function [] = allFreqBarChart(startindex, leng, js, channelnum)
%start index may have to be manually adjusted if junk trials happen
%js is the json file to be read
%channelnum is the specific channel in the lfpmontagetimedomain to be read
%color creates the color of the line

%for each frequency band we will have six bars to represent each contact
%pair
%put that with the normal plot with the clinic use one
%but only plot the coontact pair found to have max beta, just out of
%interest

ifcounter = 0;          %used to put each trial in a new column
y = zeros(4096,3);      %length fits output from pspectruum, 3 columns for each trial
for i = startindex:leng
    if strcmp(js.LfpMontageTimeDomain(i).Channel,channelnum) 
        ifcounter = ifcounter + 1;
        t = js.LfpMontageTimeDomain(i).TimeDomainData;
        [p,f] = pspectrum(t, 250, 'FrequencyLimits', [0 100]); %250 comes from json file itself
        y(:,ifcounter) = p;
    end
end 
    %label the max
    ym = mean(y, 2); %takes the mean of y by rows
    edges = [0 4; 5 8; 9 11; 12 33; 34 100];  %frequency thresholds
    ymean = zeros(5, 1);
    ystd = zeros(5, 1);
    yvar = zeros(5, 1); 
    for b = 1:size(edges, 1)
        strf = edges(b, 1);
        stpf = edges(b, 2);
        fnew = f > strf & f < stpf;
        ymean(b) = mean(ym(fnew));
        ystd(b) = std(ym(fnew));
        yvar(b) = var(ym(fnew)); 
    end
    %basically std dev values are being cwazy so! instead of using 
    %std dev, I am using var and seeing how that is better 
    figure
    bar(ymean)
    set(gca,'xticklabel',{'delta','theta','alpha','beta','gamma'})
    hold on
    errorbar(1:5, ymean, ystd) 
    titlechannel = string(replace(channelnum, '_', ' '));
    title(["(std) All frequency info for", titlechannel])
    
    figure 
    bar(ymean)
    set(gca,'xticklabel',{'delta','theta','alpha','beta','gamma'})
    hold on
    errorbar(1:5, ymean, yvar) 
    titlechannel = string(replace(channelnum, '_', ' '));
    title(["(var) All frequency info for", titlechannel])
end 