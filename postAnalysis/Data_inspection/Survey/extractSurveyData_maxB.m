function [outTABLE] = extractSurveyData_maxB(js, side, startIndex)

channels = unique({js.LfpMontageTimeDomain.Channel}, 'stable'); %should get 6 or 12 channels

%so w single battery pts, this basically does nothing bc only hemisphere
%data is present
%w double battery pts this takes into account which side to get
%also good double check that the pt directory is consistent

if contains(side, 'r',  'IgnoreCase', true) %select based on sides
    channels = channels(contains(channels, 'RIGHT'));
elseif contains(side, 'l', 'IgnoreCase', true)
    channels = channels(contains(channels, 'LEFT'));
end

leng = numel({js.LfpMontageTimeDomain.Channel});
% Time domain data
tdData = js.LfpMontageTimeDomain;
tdData = tdData(startIndex:leng);
% All channels and runs LABELS
allcAr = {tdData.Channel};
allraw = {tdData.TimeDomainData};

% 4096 samples is the number of samples generated from computing the power
% spectral info from 250Hz signal - limited to 1-100 Hz.
numChans = 6;
numRuns = 3;
numSides = 1;
outDat = cell(numChans*numRuns*numSides,1);
outRUn = nan(numChans*numRuns*numSides,1);
outSide = cell(numChans*numRuns*numSides,1);
outChan = cell(numChans*numRuns*numSides,1);

combCount = 1;

% Loop throgh Lfp Montage Time Domain
%if doubleBattery

for ci = 1:length(channels)

    chanUnique = channels{ci};
    uniCruns = allraw(matches(allcAr,chanUnique));

    % Channel string info
    chanNameParts = split(chanUnique,'_');
    chanNAME = join(chanNameParts(1:3),'_');
    %sideNAME = chanNameParts{4}(1);

    for ri = 1:length(uniCruns)

        uniRUN = uniCruns{ri};
        % Use known Fs = 250 Hz and limit to 0 - 100 Hz
        [Pxx , Fxx] = pwelch(double(uniRUN), hanning(250), 125, 256, 250, 'onesided');
        uVp_t = sqrt(Pxx).*rms(hanning(250)).*sqrt(2).*2.*250/256;

        % betaBand = pwR(freQ > 12 & freQ < 33);
        outDat{combCount,1} = table(uVp_t,Fxx,'VariableNames',{'Power','Frequency'});
        outRUn(combCount,1) = ri;
        outSide{combCount,1} = side;
        outChan{combCount,1} = chanNAME{1};

        combCount = combCount + 1;

    end
end

%else

outTABLE = table(outDat,outRUn,outSide,outChan,'VariableNames',...
    {'PF_Data','RunNum','SideID','ChanID'});
end