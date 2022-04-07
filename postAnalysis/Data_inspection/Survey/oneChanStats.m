function [oneChan] = oneChanStats(outTABLE, c) 
runs = zeros(82,3);
ids = cell(82,3);
%chop larger table to just channel of interest 
temp = outTABLE((contains(outTABLE.ChanID, c)), :);

for i = 1:3 
powerdata = temp.PF_Data{i};
%find freq vales of interest
fbeta = powerdata.Frequency > 13 & powerdata.Frequency < 33;
betapwr = powerdata.Power(fbeta);
betapwr = abs(log10(betapwr));
betapwr = downsample(betapwr, 10);
x = length(betapwr); 
id = append('Run', num2str(i)); 
runs(:,i) = betapwr; 
ids(:,i) = {id}; 
end 

run_export = ([runs(1:x), runs(x+1:2*x), runs(2*x+1:3*x)])'; 
id_export = ([ids(1:x), ids(x+1:2*x), ids(2*x+1:3*x)])'; 
Tid = array2table(id_export); 
Trun = array2table(run_export); 
T = [Tid, Trun]; 
T.Properties.VariableNames = {'Identifiers', 'Values'}; 
oneChan = T; 
end 