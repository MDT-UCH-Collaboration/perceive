% close all
% mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
% sleepPlotHelp_fun(mainLOC , num2str(9))

mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
cd(mainLOC)


%% Figure 1A - Timeline plot of single patient entire sequence
close all
% To do::
% 3. add gray block for night and day 

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
cMAP = brewermap([],"RdBu");
[reMAP] = reMapCmap(medLine,minVale,maxVale,cMAP,smSunful);
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

% Map start and stop times from EVENTS
[outIND] = findBEDinds('inbed' , evData , tmData);

xline(outIND,'-',repmat({'in bed'},length(outIND),1))


title('Patient 3')




% yline(medLine)
% Create the RED up and BLUE down with gradient effect
% https://observablehq.com/@analyzer2004/plot-gallery















































function [patDATA] = getPatDat(cID , hID , tID)

matDir = dir('*.mat');
matNames = {matDir.name};
matEls = split(matNames,'_');

cIDs1 = matEls(:,:,1);
cIDs2 = extractAfter(cIDs1,4);

hIDs = matEls(:,:,2);

tIDs1 = matEls(:,:,3);
tIDs2 = extractBefore(tIDs1,'.');

matLog = matches(cIDs2,cID) & matches(hIDs,hID) & matches(tIDs2,tID);

if matches(tID,'TimeLine')
    load(matNames{matLog},'outMAT');
    tmpDaOUT = outMAT;
else
    load(matNames{matLog},'outEVENTS');
    tmpDaOUT = outEVENTS;
end

patDATA = tmpDaOUT;

end




function [reMAP] = reMapCmap(medPoint,minVal,maxVal,cMap,datA)

% Create iMap
% Half 256 0.1 0.5 mapped to lower bound
lowBound = linspace(minVal,medPoint,128);
highBound = linspace(medPoint,maxVal,128);
iMap = [lowBound , highBound];

cMap = flipud(cMap);

reMAP = zeros(length(datA),3);
for di = 1:length(datA)

    tmpDat = datA(di);
    [~, rowID] = min(abs(iMap - tmpDat));
    tCmp = cMap(rowID,:);
    reMAP(di,:) = tCmp;

end

end



function [outIND] = findBEDinds(sTATE , evDAT , tlDAT)

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
for ei = 1:length(evDAY)

    tmpDAY = evDAY(ei);
    tmpHOUR = evHOUR(ei);
    tmpMIN = evMIN(ei);
    tldayFind = ismember(tlDAY,tmpDAY);
    tlhourFind = ismember(tlHOUR,tmpHOUR) & tldayFind;
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



