function [oneChan] = oneChanStats(outTABLE, c, flb, fub, res) 
%Calculate the stats for a single channel 

% runs = zeros(82,3);
% ids = cell(82,3);
%chop larger table to just channel of interest 
%imported value based upon what we choose lmoa 
z = fub - flb; %define length of frequency range 
z = z * 1/res; %get the length of the vector we WANT 
%length * amount of samples per second 

dummy = outTABLE.PF_Data{1}; %arbitrary data 
ffreq = dummy.Frequency > flb & dummy.Frequency < fub; %find frequency
fpwr = dummy.Power(ffreq); %find power within that frequency
dwnsamp = round(length(fpwr) / z); %vector length we have / vector length we want


temp = outTABLE((contains(outTABLE.ChanID, c)), :);

for i = 1:3 
powerdata = temp.PF_Data{i};
%find freq vales of interest
fbeta = powerdata.Frequency > flb & powerdata.Frequency < fub;
betapwr = powerdata.Power(fbeta);
betapwr = abs(log10(betapwr));
betapwr = downsample(betapwr, dwnsamp);
x = length(betapwr); 
if i == 1 
    runs = zeros(x,3);
    ids = cell(x,3);
end 

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