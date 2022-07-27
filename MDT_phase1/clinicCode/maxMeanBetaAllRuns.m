function [] = maxMeanBetaAllRuns(startindex, js, leng, channel)
% %from this function I want 3  6 x 2 tablese
% %each table will tell about the single run for each contact pair
% %the table will be ordered vertically from highest to lowest beta
% %for each run, convert to power spectrum data
% %look at that data for only the beta range
% %store the mean and std in a table
% %and display
% 
% %might have more than three runs...
% %needs to dynamically change size I think?
% %and the tables need to order themselves
% %tbl = table(Std, Max)
locallength = length(channel);
means = zeros(locallength, 3);
stds = zeros(locallength, 3);
for m=  1:locallength
    ifcounter = 0;
    for i = startindex:leng
        if strcmp(js.LfpMontageTimeDomain(i).Channel, channel{m})
            ifcounter = ifcounter +1;
            t = js.LfpMontageTimeDomain(i).TimeDomainData;
            [p,f] = pspectrum(t, 250, 'FrequencyLimits', [0 100]); %250 comes from json file itself
            fnew = f > 13 & f < 30;
            means(m, ifcounter) = mean(p(fnew));
            stds(m, ifcounter) = std(p(fnew));
        end
    end
end
Tmean = array2table(means, 'VariableNames', {'Mean R1', 'Mean R2', 'Mean R3'});
Tstd = array2table(stds, 'VariableNames', {'STD R1', 'STD R2', 'STD R3'});
Tname = cell2table(channel');
tbl = [Tname Tmean Tstd];
sortrows(tbl, {'Mean R1', 'Mean R2', 'Mean R3'}, 'descend')
end