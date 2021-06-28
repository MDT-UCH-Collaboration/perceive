function [] = graphEachTrials(startindex, js, leng, channel)


ifcounter = 0;          %used to put each trial in a new column
y = zeros(4096,3);      %length fits output from pspectruum, 3 columns for each trial
locallength = length(channel);

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
    plot(f, smy);
    title(["Three trials for", legendChan{m}])
    xlim([0 60])
    xline(13);
    xline(30);
    ifcounter = 0;
end
end

