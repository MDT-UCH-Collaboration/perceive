%given the main table (dir), add the channels that have the max power in 
%frequency set by flb (freuqnecy lower bound) and fub (frequency upper
%bound)
function [maxChanMaster] = maxChanTable(dir, flb, fub)
for i = 1:length(dir.Patient) %within a pt and hemisphere
    jsonFiles = dir.JsonSurvey{i};
    js =jsondecode(fileread(jsonFiles));
    side = dir.Hemisphere{i};
    startIndex = dir.StartIndex(i);
    num = num2str(dir.Patient(i));
    hemi = upper(dir.Hemisphere(i));
    doubleBattery = dir.DoubleBattery(i); 
    channels = unique({js.LfpMontageTimeDomain.Channel}, 'stable');
    if doubleBattery == 0
        channels = channels(contains(channels, hemi)); %select for channels w that hemisphere
    end
    leng = numel({js.LfpMontageTimeDomain.Channel});
    %run func to go into a json, go thru all the channels, total and average
    %the power, store it into a table, find the max value channel, and store
    %that channel into another table that gets appended to dir
    %sort to find max

    %generate the max frequency data
    ifcounter = 0;          %used to put each trial in a new column
    y = zeros(4096,3);      %length fits output from pspectruum, 3 columns for each trial
    locallength = length(channels);
    maxBetas = zeros(6, 1);

    for m = 1:locallength %within a single channel
        for n = startIndex:leng %within a single rune
            if strcmp(js.LfpMontageTimeDomain(n).Channel,channels{m})
                ifcounter = ifcounter + 1;
                t = js.LfpMontageTimeDomain(n).TimeDomainData;
                [p,f] = pspectrum(t, 250, 'FrequencyLimits', [0 100]); %250 comes from json file itself
                y(:,ifcounter) = p;
            end
        end
        ym = mean(y, 2); %average across runs (by rows)
        fbeta = f > flb & f < fub; %logical to only find beta region
        maxBetas(m) = max(ym(fbeta)); %find max within ym within beta %FINDS THE PEAK VALUE AND KEEPS IT
        ifcounter = 0;
    end

    %label the max and table it
    chan = rows2vars(cell2table(channels));
    maxT = [chan, array2table(maxBetas)];
    maxT = removevars(maxT, "OriginalVariableNames");
    maxT.Properties.VariableNames = {'Channel', 'Average Max Power'};

    [~, y] = sort(table2array(maxT(:,2)), 'descend');
    maxT = maxT(y,:);
    %now I know that the very first channel is the max freq channel
    maxChanTemp = {maxT.Channel{1}}; 
    %follow text protocol found in ExtractBrainSenseSurveyData
    maxChanTemp = split(maxChanTemp,'_');
    maxChanTemp = join(maxChanTemp(1:3),'_');
    if i == 1 
        maxChanMaster = maxChanTemp;
    else 
        maxChanMaster = [maxChanMaster; maxChanTemp]; 
    end
end   

end 