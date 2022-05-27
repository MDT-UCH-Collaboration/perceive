function [] = graphEachTrials(startindex, js, leng, channel, maxval)


ifcounter = 0;          %used to put each trial in a new column
y = zeros(4096,3);      %length fits output from pspectruum, 3 columns for each trial
locallength = length(channel);

%clean up contact names 
legendChan = cell(6);
oldchar = {'LEFT', 'RIGHT', '0', '1', '2', '3', '4', '5', '_', 'AND'};
newchar = {''};
oldnum = {'ZERO', 'ONE', 'TWO', 'THREE'};
newnum = {'0', '1', '2', '3'};
for b = 1:length(channel)
    legendChan{b} = replace(channel{b}, oldchar, newchar);
    legendChan{b} = replace(legendChan{b}, oldnum, newnum);
end

for m = 1:locallength
    subplot(3, 3, m);
    for g = startindex:leng
        if strcmp(js.LfpMontageTimeDomain(g).Channel,channel{m})
            ifcounter = ifcounter + 1;
            t = js.LfpMontageTimeDomain(g).TimeDomainData;
            [p,f] = pspectrum(t, 250, 'FrequencyLimits', [0 100]); %250 comes from json file itself
            y(:,ifcounter) = p;
        end
    end
    smy = smoothdata(y, 1, "sgolay", 250); %y is now smooth
    %max = max(y, [], 1);
    plot(f, smy);
    title([legendChan{m}])
    xlim([0 40])
    ylim([0 maxval*1.1])
    xline(13);
    xline(30);
    ifcounter = 0;
    sgtitle("Run Comparison by Channel")
    xticks([0 20 40])
    yticks([0 round(maxval*0.5, 2) round(maxval,2)])
    xlabel("Frequency (Hz)")
    ylabel("uVp")
    
    
end
end

