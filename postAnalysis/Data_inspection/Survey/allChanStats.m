function [allChan] = allChanStats(outTABLE) 
%take each channel 
%get all the runs data 
%average it 
%put the average into a column 
%label that data with a letter 
chan = unique(outTABLE.ChanID, 'stable'); 
runs = zeros(82,3);
ids = cell(82,1);
counter = 0; 
letter = {'A', 'B', 'C', 'D', 'E', 'F'};
%chop larger table to just channel of interest 
for k = 1:length(chan)
temp = outTABLE(contains(outTABLE.ChanID, chan{k}) & contains(outTABLE.ChanID, chan{k}), :);
counter = counter +1; 

for i = 1:3 
powerdata = temp.PF_Data{i};
%find freq vales of interest
fbeta = powerdata.Frequency > 13 & powerdata.Frequency < 33;
betapwr = powerdata.Power(fbeta);
betapwr = abs(log10(betapwr));
betapwr = downsample(betapwr, 10);
x = length(betapwr);  
runs(:,i) = betapwr; 
end 

run_export = mean(runs, 2); 
ids = repmat({letter(counter)}, 82, 1);
Tid = array2table(ids); 
Trun = array2table(run_export); 
T = [Tid, Trun]; 
if counter == 1 
    Tmaster = T; 
else 
    Tmaster = [Tmaster; T]; 
end  
end
Tmaster.Properties.VariableNames = {'Identifiers', 'Values'}; 
allChan = Tmaster;
end 