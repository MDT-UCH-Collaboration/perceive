function [allChan] = allChanStats(outTABLE, flb, fub, res)
%take each channel
%get all the runs data
%average it
%put the average into a column
%label that data with a letter
z = fub - flb; %define length of frequency range
z = z * 1/res; %get the length of the vector we WANT
%length * amount of samples per second

dummy = outTABLE.PF_Data{1}; %arbitrary data
ffreq = dummy.Frequency > flb & dummy.Frequency < fub; %find frequency
fpwr = dummy.Power(ffreq); %find power within that frequency
dwnsamp = round(length(fpwr) / z); %vector length we have / vector length we want

chan = unique(outTABLE.ChanID, 'stable');
runs = zeros(82,3);
ids = cell(82,1);
counter = 0;
letter = {'A', 'B', 'C', 'D', 'E', 'F'};
%chop larger table to just channel of interest
datatab = [];
for k = 1:length(chan)
    temp = outTABLE(contains(outTABLE.ChanID, chan{k}) & contains(outTABLE.ChanID, chan{k}), :);
    counter = counter +1;

    for i = 1:3
        powerdata = temp.PF_Data{i};
        %find freq vales of interest
        fbeta = powerdata.Frequency > flb & powerdata.Frequency < fub;
        betapwr = powerdata.Power(fbeta);
%         betapwr = abs(log10(betapwr)); % don't log 10
        betapwr = downsample(betapwr, dwnsamp);
        x = length(betapwr);
        runtab = repmat(temp.RunNum(i), x, 1);
        chantab = repmat(temp.ChanID(i), x, 1);
        localtab = array2table(betapwr);
        localtab.Run = runtab;
        localtab.Chan = chantab;
        datatab = [datatab; localtab];
    end

end
allChan = datatab;
end