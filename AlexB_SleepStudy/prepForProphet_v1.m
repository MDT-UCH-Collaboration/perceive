function [] = prepForProphet_v1()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
cd(mainLOC)
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy');
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy\matplotlib03232022');

subjectID = '3';
hemisphere = 'L';
[tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');
% [apData] = getPatDat(subjectID , hemisphere , 'ActALL');
% [cleanApT] = trim144Apdata(tmData , apData);

allTIME = tmData.timeDATA;

% Clean up time data



nLFP = tmData.LFP;
unfurlLFP = nLFP(:);
mSunful = unfurlLFP - (min(unfurlLFP));
mSunful(mSunful > 2.3000e+09) = nan;
mSunful = normalize(mSunful, 'range');
smSunful = smoothdata(mSunful,'rloess',20,'omitnan');


end





function [patDATA] = getPatDat(cID , hID , tID)

matDir = dir('*.mat');
matNames = {matDir.name};
matEls = split(matNames,'_');

cIDs1 = matEls(:,:,1);
cIDs2 = extractAfter(cIDs1,4);

hIDs = matEls(:,:,2);

tIDs1 = matEls(:,:,3);
tIDs2 = extractBefore(tIDs1,'.');

if matches(tID,'TimeLine')
    matLog = matches(cIDs2,cID) & matches(hIDs,hID) & matches(tIDs2,tID);
    load(matNames{matLog},'outMAT');
    tmpDaOUT = outMAT;
elseif matches(tID,'Events')
    matLog = matches(cIDs2,cID) & matches(hIDs,hID) & matches(tIDs2,tID);
    load(matNames{matLog},'outEVENTS');
    tmpDaOUT = outEVENTS;
else
    matLog = matches(cIDs2,cID) & matches(tIDs2,tID);
    load(matNames{matLog},'rawActSlWk')
    tmpDaOUT = rawActSlWk;
end

patDATA = tmpDaOUT;

end




function [cleanApT] = trim144Apdata(proCombdat , actPydat)

allMonth = month(datetime(join([actPydat.Date actPydat.Time],' ')));
allDay = day(datetime(join([actPydat.Date actPydat.Time],' ')));
allHour = hour(datetime(join([actPydat.Date actPydat.Time],' ')));
allMin = minute(datetime(join([actPydat.Date actPydat.Time],' ')));
allSec = second(datetime(join([actPydat.Date actPydat.Time],' ')));

% loop through row and column
cosinar = nan(size(proCombdat.actTime));
fractaL = nan(size(proCombdat.actTime));
SSA = nan(size(proCombdat.actTime));
activity = nan(size(proCombdat.actTime));
cirCad = nan(size(proCombdat.actTime));
ultrid = nan(size(proCombdat.actTime));
recSt = nan(size(proCombdat.actTime));
ckSW = nan(size(proCombdat.actTime)); % 0: sleep, 1: wake
sadehSW = nan(size(proCombdat.actTime));
scripsSW = nan(size(proCombdat.actTime));
oakleySW = nan(size(proCombdat.actTime));
crespoSW = nan(size(proCombdat.actTime));
ronenbSW = nan(size(proCombdat.actTime));

for dayI = 1:size(proCombdat.actTime,2)

    tmpMonth1 = unique(month(proCombdat.actTime(:,dayI)));
    tmpMonth2 = tmpMonth1(~isnan(tmpMonth1));
    tmpDay1 = unique(day(proCombdat.actTime(:,dayI)));
    tmpDay2 = tmpDay1(~isnan(tmpDay1));

    mdCheck = allMonth == tmpMonth2 & allDay == tmpDay2;
    if sum(mdCheck) == 0 % if day is missing
        % Fill column with nan, so ski;
        continue
    else
        for hourMin = 1:size(proCombdat.actTime,1)
            tmpHou1 = unique(hour(proCombdat.actTime(hourMin,dayI)));
            tmpMin1 = unique(minute(proCombdat.actTime(hourMin,dayI)));
            tmpSec = unique(second(proCombdat.actTime(hourMin,dayI)));
            if tmpSec > 30
                tmpSec = 30;
            else
                tmpSec = 0;
            end
            % If on last row contine
            chkMIN = find((tmpHou1 == allHour &...
                tmpMin1 == allMin & tmpSec == allSec & mdCheck));
            if sum(chkMIN) == 0 || isempty(chkMIN)
                continue
            else

                % Fits

                cosinar(hourMin,dayI) = actPydat.Cosine(chkMIN);
                fractaL(hourMin,dayI) = actPydat.Fractal(chkMIN);
                SSA(hourMin,dayI)     = actPydat.SSA(chkMIN);
                cirCad(hourMin,dayI)  = actPydat.Circadian(chkMIN);
                ultrid(hourMin,dayI)  = actPydat.Ultridian(chkMIN);
                recSt(hourMin,dayI)   = actPydat.Reconstruct(chkMIN);
                activity(hourMin,dayI) = str2double(actPydat.Activity{chkMIN});
                % SW
                ckSW(hourMin,dayI) = actPydat.CK(chkMIN);
                sadehSW(hourMin,dayI) = actPydat.Sadeh(chkMIN);
                scripsSW(hourMin,dayI) = actPydat.Scripps(chkMIN);
                oakleySW(hourMin,dayI) = actPydat.Oakley(chkMIN);
                crespoSW(hourMin,dayI) = actPydat.Crespo(chkMIN);
                ronenbSW(hourMin,dayI) = actPydat.Roneenberg(chkMIN);

            end % If/else check if there are minutes/hours

        end % Loop through 10 min blocks
    end % If/else check if there are days
end % Loop through days

% Create output
cleanApT.cosinar = cosinar;
cleanApT.fractaL = fractaL;
cleanApT.SSA = SSA;
cleanApT.cirCad = cirCad;
cleanApT.ultrid = ultrid;
cleanApT.recSt = recSt;
cleanApT.Activity = activity;
cleanApT.ckSW = ckSW;
cleanApT.sadehSW = sadehSW;
cleanApT.scripsSW = scripsSW;
cleanApT.oakleySW = oakleySW;
cleanApT.crespoSW = crespoSW;
cleanApT.ronenbSW = ronenbSW;


end
