function [] = createAlexPerceptFigures(figNUM)

% mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
% sleepPlotHelp_fun(mainLOC , num2str(9))

mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
cd(mainLOC)
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy');
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy\matplotlib03232022');

%% Figure 1A - Timeline plot of single patient entire sequence - LFP
close all
% To do::
% 3. Scale color to day/night cycle and not the median

switch figNUM
    case 1

        subjectID = '3';
        hemisphere = 'L';
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');
        [evData] = getPatDat(subjectID , hemisphere , 'Events');

        fig = figure;
        left_color = [1 0 0];
        right_color = [0 0 0];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);

        nLFP = tmData.LFP;
        unfurlLFP = nLFP(:);
        mSunful = unfurlLFP - (min(unfurlLFP));
        mSunful(mSunful > 2.3000e+09) = nan;
        mSunful = normalize(mSunful, 'range');
        smSunful = smoothdata(mSunful,'rloess',20,'omitnan');

        % p2 = scatter(1:length(mSunful),mSunful,'filled','o','ColorVariable','r','SizeData',20);
        % alpha(p2 , 0.3)
        medLine = median(smSunful);
        minVale = min(smSunful);
        maxVale = max(smSunful);
        % iMAP = linspace(0,1,256);
        % cMAP = brewermap([],"RdBu");
        cMAP = cividis;
        [reMAP] = reMapCmap(medLine,minVale,maxVale,cMAP,smSunful,0);
        scatter(1:length(smSunful),smSunful,[],reMAP,'filled')
        ylim([0 round(maxVale + 0.1,1)])
        yticks([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
        yticklabels([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
        ylabel('Scaled power')
        dayStarts = round(linspace(1,length(smSunful)-144,length(smSunful)/144));
        xticks(dayStarts);
        xticklabels(1:length(smSunful)/144)
        xlim([1, length(smSunful)])
        xlabel('Days of recording')
        set(gca,'TickLength',[0 .001])

        %%%% GET IN BED times Map start and stop times from EVENTS
        [inBEDx , inBEDdays] = findBEDinds('inbed' , evData , tmData);

        % Need line code
        % For IN BED
        % 1. yMIN = yValue
        % 2. yMAX = max Y axis
        hold on
        yLIMn = ylim();
        inBEDlineX = transpose([inBEDx , inBEDx]);
        inBEDminY = smSunful(inBEDx) + smSunful(inBEDx)*0.1;
        inBEDmaxY = repmat(yLIMn(2),length(inBEDx),1);
        inBEDlineY = transpose([inBEDminY ,  inBEDmaxY]);
        line(inBEDlineX , inBEDlineY , 'Color', 'k', 'LineWidth',1);

        % IN Bed text code - Right justified
        xt_IN = inBEDx + 2;
        yt_IN = repmat(yLIMn(2)*0.95,length(inBEDx),1);
        str_IN = repmat({'IN'},length(inBEDx),1);
        inText = text(xt_IN,yt_IN,str_IN);
        % set(inText,'Rotation',90)
        % xline(inBEDx,'-',repmat({'in bed'},length(inBEDx),1))

        %%%% GET IN BED times Map start and stop times from EVENTS
        [outBEDx , outBEDdays] = findBEDinds('wake' , evData , tmData);

        % Out BED line code
        % For IN BED
        % 1. yMIN = yValue
        % 2. yMAX = max Y axis
        outBEDlineX = transpose([outBEDx , outBEDx]);
        outBEDminY = smSunful(outBEDx) - smSunful(outBEDx)*0.1;
        outBEDmaxY = repmat(yLIMn(1),length(outBEDx),1);
        outBEDlineY = transpose([outBEDminY ,  outBEDmaxY]);
        line(outBEDlineX , outBEDlineY , 'Color', 'k', 'LineWidth',1);

        % OUT Bed text code - Right justified
        xt_OUT = outBEDx - 2;
        yt_OUT = repmat(yLIMn(2) - (yLIMn(2)*0.95),length(outBEDx),1);
        str_OUT = repmat({'OUT'},length(outBEDx),1);
        outText = text(xt_OUT,yt_OUT,str_OUT,'HorizontalAlignment','right');

        % Add grey patch code [above and below]
        % First find matching IN and OUT of bed days
        % ASSUMING IN BED prior to MIDNIGHT AND OUT OF BED AFTER MIDNIGHT next DAY
        inANDout = nan(length(inBEDdays),2);
        inOUTc = 1;
        for ib = 1:length(inBEDdays)

            tmpINbd = inBEDdays(ib);
            tmpOUTdb = tmpINbd + 1;
            outLog = ismember(outBEDdays,tmpOUTdb);

            if sum(outLog)
                inANDout(inOUTc,1) = inBEDx(ib);
                inANDout(inOUTc,2) = outBEDx(outLog);
                inOUTc = inOUTc + 1;
            else
                if tmpINbd == 30 || tmpINbd == 31
                    tmpOUTdb = 1;
                    outLog = ismember(outBEDdays,tmpOUTdb);
                    inANDout(inOUTc,1) = inBEDx(ib);
                    inANDout(inOUTc,2) = outBEDx(outLog);
                    inOUTc = inOUTc + 1;
                end
            end
        end

        inANDouti = inANDout(~isnan(inANDout(:,1)),:);

        % ADD patch
        for ip = 1:size(inANDouti,1)

            xPbot = inANDouti(ip,1):inANDouti(ip,2);
            xPtop = fliplr(xPbot);
            yPbot = repmat(yLIMn(1),size(xPbot));
            yPtop = fliplr(transpose(smSunful(xPbot)));

            patch([xPbot xPtop],[yPbot yPtop],[0.5 0.5 0.5],...
                'FaceAlpha',0.4,'EdgeColor','none')




        end

        title('Patient 3')

    case 2 % Figure 1B - Timeline plot of single patient entire sequence - ACTIGRAPHY
        % To do::
        % 1. 
        % 2. 
        % 3. Sleep Wake
        % 4. Plot raw
        % 5. Plot fit

        subjectID = '3';
        hemisphere = 'L';
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');
        [apData] = getPatDat(subjectID , hemisphere , 'ActALL');
        [cleanApT] = trim144Apdata(tmData , apData);

        test = 1;

        figure;
        yyaxis left
        plot(cleanApT.cirCad(:))
        hold on
        yyaxis right
        plot(tmData.ActMean(:))
%         ylim([0 1.5])

        yyaxis left
        plot(apData.Circadian)
        hold on
        yyaxis right
        plot(str2double(apData.Activity))



        % Figure 1C - Heat map of actigraphy/lfp






        % Figure 1D - Event plot for rep patient



























end

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




function [reMAP] = reMapCmap(medPoint,minVal,maxVal,cMap,datA,flipFlag)

% Create iMap
% Half 256 0.1 0.5 mapped to lower bound
lowBound = linspace(minVal,medPoint,128);
highBound = linspace(medPoint,maxVal,128);
iMap = [lowBound , highBound];

if flipFlag
    cMap = flipud(cMap);
end

reMAP = zeros(length(datA),3);
for di = 1:length(datA)

    tmpDat = datA(di);
    [~, rowID] = min(abs(iMap - tmpDat));
    tCmp = cMap(rowID,:);
    reMAP(di,:) = tCmp;

end

end



function [outIND , dayOUTall] = findBEDinds(sTATE , evDAT , tlDAT)

if matches(sTATE,'inbed')
    evDAY = evDAT.GoingToBed.Day;
    evHOUR = evDAT.GoingToBed.Hour;
    evMIN = evDAT.GoingToBed.Minute;
else
    evDAY = evDAT.WakingUp.Day;
    evHOUR = evDAT.WakingUp.Hour;
    evMIN = evDAT.WakingUp.Minute;
end

tlDAY = tlDAT.day(:);
tlHOUR = tlDAT.hour(:);
tlMIN = tlDAT.minu(:);

outIND = zeros(length(evDAY),1);
dayOUTall = zeros(length(evDAY),1);
for ei = 1:length(evDAY)

    tmpDAY = evDAY(ei);
    tmpHOUR = evHOUR(ei);
    tmpMIN = evMIN(ei);
    tldayFind = ismember(tlDAY,tmpDAY);
    tlhourFind = ismember(tlHOUR,tmpHOUR) & tldayFind;
    dayOUTall(ei) = tmpDAY;
    tlminFind = find(ismember(tlMIN,tmpMIN) & tlhourFind);
    if isempty(tlminFind)
        hourINDS = find(tlhourFind);
        hourTRIM = tlMIN(tlhourFind);
        [~ , minLOC] = min(abs(hourTRIM - tmpMIN));
        tlminFind = hourINDS(minLOC);
    end
    outIND(ei) = tlminFind;

end


end



function [cleanApT] = trim144Apdata(proCombdat , actPydat)

allMonth = month(datetime(join([actPydat.Date actPydat.Time],' ')));
allDay = day(datetime(join([actPydat.Date actPydat.Time],' ')));
allHour = hour(datetime(join([actPydat.Date actPydat.Time],' ')));
allMin = minute(datetime(join([actPydat.Date actPydat.Time],' ')));
% allSec = second(datetime(join([actPydat.Date actPydat.Time],' ')));

% loop through row and column
cosinar = nan(size(proCombdat.actTime));
fractaL = nan(size(proCombdat.actTime));
SSA = nan(size(proCombdat.actTime));
cirCad = nan(size(proCombdat.actTime));
ultrid = nan(size(proCombdat.actTime));
recSt = nan(size(proCombdat.actTime));
ckSW = nan([size(proCombdat.actTime),2]); % 0: sleep, 1: wake
sadehSW = nan([size(proCombdat.actTime),2]);
scripsSW = nan([size(proCombdat.actTime),2]);
oakleySW = nan([size(proCombdat.actTime),2]);
crespoSW = nan([size(proCombdat.actTime),2]);
ronenbSW = nan([size(proCombdat.actTime),2]); 

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
                
                % tmpSec = unique(second(proCombdat.actTime(hourMin,dayI)));
                % If on last row contine

                chkMIN = find((tmpHou1 == allHour &...
                    tmpMin1 == allMin & mdCheck));
                if sum(chkMIN) == 0
                    continue
                else

                    startMIN = find((tmpHou1 == allHour &...
                        tmpMin1 == allMin & mdCheck),1,'first');
                    if hourMin == height(proCombdat.actTime)
                        stopMIN = find((tmpHou1 == allHour &...
                            tmpMin1 == allMin & mdCheck),1,'last');
                    else
                        tmpHou2 = unique(hour(proCombdat.actTime(hourMin + 1,dayI)));
                        tmpMin2 = unique(minute(proCombdat.actTime(hourMin + 1,dayI)));
                        stopMIN = find((tmpHou2 == allHour &...
                            tmpMin2 == allMin & mdCheck),1,'first');
                    end

                    % Fits
                    cosinar(hourMin,dayI) = mean(actPydat.Cosine(startMIN:stopMIN),'omitnan');
                    fractaL(hourMin,dayI) = mean(actPydat.Fractal(startMIN:stopMIN),'omitnan');
                    SSA(hourMin,dayI)     = mean(actPydat.SSA(startMIN:stopMIN),'omitnan');
                    cirCad(hourMin,dayI)  = mean(actPydat.Circadian(startMIN:stopMIN),'omitnan');
                    ultrid(hourMin,dayI)  = mean(actPydat.Ultridian(startMIN:stopMIN),'omitnan');
                    recSt(hourMin,dayI)   = mean(actPydat.Reconstruct(startMIN:stopMIN),'omitnan');
                    % SW
                    [ckfrac, ckState] = getSWstate(actPydat.CK,startMIN,stopMIN);
                    ckSW(hourMin,dayI,1) = ckfrac;
                    ckSW(hourMin,dayI,2) = ckState;

                    [sadfrac, sadState] = getSWstate(actPydat.Sadeh,startMIN,stopMIN);
                    sadehSW(hourMin,dayI,1) = sadfrac;
                    sadehSW(hourMin,dayI,2) = sadState;

                    [scifrac, sciState] = getSWstate(actPydat.Scripps,startMIN,stopMIN);
                    scripsSW(hourMin,dayI,1) = scifrac;
                    scripsSW(hourMin,dayI,2) = sciState;

                    [oakfrac, oakState] = getSWstate(actPydat.Oakley,startMIN,stopMIN);
                    oakleySW(hourMin,dayI,1) = oakfrac;
                    oakleySW(hourMin,dayI,2) = oakState;

                    [crefrac, creState] = getSWstate(actPydat.Crespo,startMIN,stopMIN);
                    crespoSW(hourMin,dayI,1) = crefrac;
                    crespoSW(hourMin,dayI,2) = creState;

                    [roefrac, roeState] = getSWstate(actPydat.Roneenberg,startMIN,stopMIN);
                    ronenbSW(hourMin,dayI,1) = roefrac;
                    ronenbSW(hourMin,dayI,2) = roeState;

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
cleanApT.ckSW_fracW = ckSW(:,:,1);
cleanApT.ckSW_st = ckSW(:,:,2);
cleanApT.sadehSW_fracW = sadehSW(:,:,1);
cleanApT.sadehSW_st = sadehSW(:,:,2);
cleanApT.scripsSW_fracW = scripsSW(:,:,1);
cleanApT.scripsSW_st = scripsSW(:,:,2);
cleanApT.oakleySW_fracW = oakleySW(:,:,1);
cleanApT.oakleySW_st = oakleySW(:,:,2);
cleanApT.crespoSW_fracW = crespoSW(:,:,1);
cleanApT.crespoSW_st = crespoSW(:,:,2);
cleanApT.ronenbSW_fracW = ronenbSW(:,:,1);
cleanApT.ronenbSW_st = ronenbSW(:,:,2);

end



function [outFrac, outState] = getSWstate(inSTATE,start,stop)

ckfrac = sum(inSTATE(start:stop))/length(inSTATE(start:stop));
outFrac = ckfrac;
if outFrac <= 0.4
    outState = 0;
else
    outState = 1;
end

end
