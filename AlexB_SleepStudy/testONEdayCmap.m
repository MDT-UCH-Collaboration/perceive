function [] = testONEdayCmap()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
cd(mainLOC)
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy');
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy\matplotlib03232022');
close all

% Set up Figure window
% 4 ROWS | 4 COLUMNS
% ROW 1 - Cols 1-3 = LFP plot | Col 4 Heat plot
mainFig = figure;
set(mainFig,'Position', [1485 311 1128 863]);
tiledlayout(1,4,"Padding","tight");
% Load data
subjectID = '3';
hemisphere = 'L';
[tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');
[evData] = getPatDat(subjectID , hemisphere , 'Events');

nLFP = tmData.LFP;
unfurlLFP = nLFP(:);
mSunful = unfurlLFP - (min(unfurlLFP));
mSunful(mSunful > 2.3000e+09) = nan;
mSunful = normalize(mSunful, 'range');
smSunful = smoothdata(mSunful,'rloess',20,'omitnan');
firstDayTime6At6p = [find(tmData.hour(:,2) == 6,1,'first')+144 ...
    find(tmData.hour(:,2) == 18,1,'first')+144] ;

maxVale = max(smSunful);

% Plot 1 ##########################################################
nexttile([1 4])

cMAP = cividis;
[reMAP] = reMapCmap(tmData,cMAP,smSunful,1,'timeBased');
lightCM = reMAP(10,:);
darkCM = reMAP(246,:);
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

[inBEDx , inBEDdays] = findBEDinds('inbed' , evData , tmData);

hold on
yLIMn = ylim();
inBEDlineX = transpose([inBEDx , inBEDx]);
inBEDminY = smSunful(inBEDx) + smSunful(inBEDx)*0.1;
inBEDmaxY = repmat(yLIMn(2),length(inBEDx),1);
inBEDlineY = transpose([inBEDminY ,  inBEDmaxY]);
line(inBEDlineX , inBEDlineY , 'Color', 'k', 'LineWidth',1);

% Add 'IN' to only first line
xt_IN = inBEDx(1) - 35;
yt_IN = yLIMn(2)*0.95;
str_IN = 'IN';
inText = text(xt_IN,yt_IN,str_IN);
inText.FontWeight = "bold";

[outBEDx , outBEDdays] = findBEDinds('wake' , evData , tmData);

outBEDlineX = transpose([outBEDx , outBEDx]);
outBEDminY = smSunful(outBEDx) - smSunful(outBEDx)*0.1;
outBEDmaxY = repmat(yLIMn(1),length(outBEDx),1);
outBEDlineY = transpose([outBEDminY ,  outBEDmaxY]);
line(outBEDlineX , outBEDlineY , 'Color', 'k', 'LineStyle','-.','LineWidth',1);

xt_OUT = outBEDx(1) + 5;
yt_OUT = yLIMn(2) - (yLIMn(2)*0.95);
str_OUT = 'OUT';
outText = text(xt_OUT,yt_OUT,str_OUT);
outText.FontWeight = "bold";

% Add 6AM to 6PM markers
xl6a = xline(firstDayTime6At6p(1),'--','6AM','LabelVerticalAlignment','bottom',...
    'LabelHorizontalAlignment','center');
xl6a.Color = cMAP(246,:);
xl6p = xline(firstDayTime6At6p(2),'--','6PM','LabelVerticalAlignment','bottom',...
    'LabelHorizontalAlignment','center');
xl6p.Color = cMAP(10,:);

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





function [reMAP] = reMapCmap(inDATA,cMap,datA,flipFlag,reMapS)

% Create iMap
% Half 256 0.1 0.5 mapped to lower bound

switch reMapS
    case 'median'

        medPoint = median(inDATA);
        minVal = min(inDATA);
        maxVal = max(inDATA);

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

    case 'ActSleepWake'



    case 'timeBased'

        timMat1 = repmat(linspace(0,0.5,6),24,1) + transpose(0:23);
        timMat2 = timMat1(:);
        timMat3 = sort(timMat2);
        INBEDave = 8; % CHECK per patient
        OUTBEDave = 20;
        timMatf = [timMat3(timMat3 >= INBEDave & timMat3 <= OUTBEDave) ; ...
            timMat3(timMat3 > OUTBEDave) ; timMat3(timMat3 < INBEDave)];

        minMat = inDATA.minu;
        houMat = inDATA.hour;
        % Consistify minute data
        for mi = 1:size(minMat,1)
            tmpRm = minMat(mi,:);
            uRm = unique(tmpRm);

            tmpRh = houMat(mi,:);
            uRh = unique(tmpRh);

            maxCm = 0;
            maxIDm = 0;
            for um = 1:length(uRm)
                rowS = sum(ismember(tmpRm,uRm(um)));
                if rowS > maxCm
                    maxCm = rowS;
                    maxIDm = uRm(um);
                end
            end
            minMat(mi,:) = maxIDm;

            maxCh = 0;
            maxIDh = 0;
            for uh = 1:length(uRh)
                rowS = sum(ismember(tmpRh,uRh(uh)));
                if rowS > maxCh
                    maxCh = rowS;
                    maxIDh = uRh(uh);
                end
            end
            houMat(mi,:) = maxIDh;
        end

        % unfurl hour and min and line up
        minMatu = minMat(:);
        houMatu = houMat(:);

        combMHu = (minMatu/100) + houMatu;

        if flipFlag
            cMap = flipud(cMap);
        end

        reMAP = zeros(length(datA),3);
        for di = 1:length(datA)

            tmpDat = combMHu(di);

            tmLoc = find(timMatf == tmpDat);

            cmapLoc = round((tmLoc/144)*256);

            tCmp = cMap(cmapLoc,:);
            reMAP(di,:) = tCmp;

        end

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