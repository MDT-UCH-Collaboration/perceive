function [] = makeAlexPerceptFigures_v1(figureNUM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
cd(mainLOC)
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy');
addpath('C:\Users\John\Documents\GitHub\perceive\AlexB_SleepStudy\matplotlib03232022');
close all

% To do
%. 1 Fix Figure 1A
%. 2 Add x-label to 1E [Frequency Hz]
%. 3 Add xline d t a b g

switch figureNUM
    case 1

        % Set up Figure window
        % 4 ROWS | 4 COLUMNS
        % ROW 1 - Cols 1-3 = LFP plot | Col 4 Heat plot
        mainFig = figure;
        set(mainFig,'Position', [1485 311 1128 863]);
        tiledlayout(3,4,"Padding","tight");
        % Load data
        subjectID = '3';
        hemisphere = 'L';
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');
        [evData] = getPatDat(subjectID , hemisphere , 'Events');
        [apData] = getPatDat(subjectID , hemisphere , 'ActALL');
        [cleanApT] = trim144Apdata(tmData , apData);

        nLFP = tmData.LFP;
        unfurlLFP = nLFP(:);
        mSunful = unfurlLFP - (min(unfurlLFP));
        mSunful(mSunful > 2.3000e+09) = nan;
        mSunful = normalize(mSunful, 'range');
        smSunful = smoothdata(mSunful,'rloess',20,'omitnan');

        medLine = median(smSunful);
        minVale = min(smSunful);
        maxVale = max(smSunful);

        % Plot 1 ##########################################################
        nexttile([1 4])

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

        [inBEDx , inBEDdays] = findBEDinds('inbed' , evData , tmData);

        hold on
        yLIMn = ylim();
        inBEDlineX = transpose([inBEDx , inBEDx]);
        inBEDminY = smSunful(inBEDx) + smSunful(inBEDx)*0.1;
        inBEDmaxY = repmat(yLIMn(2),length(inBEDx),1);
        inBEDlineY = transpose([inBEDminY ,  inBEDmaxY]);
        line(inBEDlineX , inBEDlineY , 'Color', 'k', 'LineWidth',1);

        xt_IN = inBEDx + 2;
        yt_IN = repmat(yLIMn(2)*0.95,length(inBEDx),1);
        str_IN = repmat({'IN'},length(inBEDx),1);
        inText = text(xt_IN,yt_IN,str_IN);

        [outBEDx , outBEDdays] = findBEDinds('wake' , evData , tmData);

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

        % PLOT 2 ##########################################################
        nexttile([1 4])

        subjectID = '3';
        hemisphere = 'L';
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');
        [apData] = getPatDat(subjectID , hemisphere , 'ActALL');
        [cleanApT] = trim144Apdata(tmData , apData);
        [evData] = getPatDat(subjectID , hemisphere , 'Events');

        % Activity
        nActraw = cleanApT.Activity;
        unfurlACTr = nActraw(:);
        sMnRmACTr = smoothdata(unfurlACTr,'gaussian',40,'omitnan');
        nRmACTr = normalize(sMnRmACTr, 'range');
        lraw = plot(nRmACTr,'LineWidth',2);
        lraw.Color = 'k';
        lraw.LineStyle = "-";

        hold on
        % Circadian Fit
        nFitraw = cleanApT.cirCad;
        unfurlFITr = nFitraw(:);
        sMnRmFITr = smoothdata(unfurlFITr,'gaussian',40,'omitnan');
        nRmFITr = normalize(sMnRmFITr, 'range');
        lcosin = plot(nRmFITr,'LineWidth',2);
        lcosin.Color = [0.5 0.5 0.5];
        lcosin.LineStyle = "-.";

        % Sleep Wake State
        % Ronnenberg: 1 = Sleep
        % Crespo: 1 = Wake
        % Invert Ronneberg
        reonUF = cleanApT.ronenbSW(:);
        nonNanInd1 = ~isnan(reonUF);
        reonUnfurli = ~reonUF(nonNanInd1);
        reonUFi = reonUF;
        reonUFi(nonNanInd1) = reonUnfurli;
        % Get Crespo
        cresUF = cleanApT.crespoSW(:);
        % Find agreement
        pairRC = [reonUFi , cresUF];
        nanInd2 = isnan(pairRC(:,1));
        % Find nonNans
        pairMatch = pairRC(:,1) == pairRC(:,2);
        swFinMat = pairRC(:,1);
        swFinMat(pairMatch) = pairRC(pairMatch,1);
        swFinMat(~pairMatch) = nan;
        swFinMat(nanInd2) = nan;

        xAxisAct = 1:length(swFinMat);
        p1 = plot(xAxisAct(swFinMat == 1),ones(size(xAxisAct(swFinMat == 1))));
        p1.LineStyle = "none";
        p1.Color = cMAP(246,:);
        p1.Marker = 'o';
        p2 = plot(xAxisAct(swFinMat == 0),zeros(size(xAxisAct(swFinMat == 0))));
        p2.Color = cMAP(10,:);
        p2.Marker = 'o';
        p2.LineStyle = "none";
        ylim([-0.1 1.1])
        leg2 = legend('Raw activity','Cosinor','Wake','Sleep');
        leg2.Position = [0.0490 0.6406 0.1126 0.0713];

        dayStarts = round(linspace(1,length(swFinMat)-144,length(swFinMat)/144));
        xticks(dayStarts);
        xticklabels(1:length(swFinMat)/144)
        xlim([1, length(swFinMat)])
        xlabel('Days of recording')
        set(gca,'TickLength',[0 .001])
        ylabel('Scaled activity')
        yticks(linspace(0,1,3))

        % Plot 3 ##########################################################
        % Cross Correlation
        nexttile(9)
        nActraw = cleanApT.Activity;
        unfurlACTr = nActraw(:);
        sMnRmACTr = smoothdata(unfurlACTr,'gaussian',40,'omitnan');
        nRmACTr = normalize(sMnRmACTr, 'range');

        nonNanLocs = ~isnan(nRmACTr);

        [crossCorvals,lags] = xcorr(smSunful(nonNanLocs),nRmACTr(nonNanLocs));
        crossCorvalsN = crossCorvals/max(crossCorvals);
        [~,ccnLoc] = max(crossCorvalsN);
        maxClagLoc = lags(ccnLoc);

        plot(lags,crossCorvalsN,'Color','k','LineWidth',1.5)
        xl1 = xline(maxClagLoc,'-.',['Max Lag = ' num2str(maxClagLoc)]);
        xl1.LineWidth = 1.5;
        xl1.Color = 'r';
        xl1.LabelVerticalAlignment = "bottom";
%         text(maxClagLoc + 100, 0.5 ,['Lag: ' int2str(maxClagLoc)])
        ylabel('\textit{r}','Interpreter','latex')
        yticks(linspace(0,1,3))
        xlabel('Lag in samples')
        title('Max lag between activity and LFP')
        axis tight

        % Plot 4 ##########################################################
        % Day v Night cross corr estimates
        nexttile(10)

%         [R,P,~,~] = corrcoef(smSunful(nonNanLocs),nRmACTr(nonNanLocs));
        % Day vs night
        % Invert Ronneberg
        reonUF = cleanApT.ronenbSW(:);
        nonNanInd1 = ~isnan(reonUF);
        reonUnfurli = ~reonUF(nonNanInd1);
        reonUFi = reonUF;
        reonUFi(nonNanInd1) = reonUnfurli;
        % Get Crespo
        cresUF = cleanApT.crespoSW(:);
        % Find agreement
        pairRC = [reonUFi , cresUF];
        % Find nonNans
        nanInd2 = isnan(pairRC(:,1));
        %         pairNnans = pairRC(nonNanInd2,:);
        pairMatch = pairRC(:,1) == pairRC(:,2);
        swFinMat = pairRC(:,1);
        swFinMat(pairMatch) = pairRC(pairMatch,1);
        swFinMat(~pairMatch) = nan;
        swFinMat(nanInd2) = nan;
        [dayBlocks , nightBlocks] = getBlocks(swFinMat);
        % 1. Pull out night and day points
        % 2. Loop through
        % 3. Correlation coefficient for each R/P
        % 4. Store by segment
        nightData = nan(length(nightBlocks),2);
        dayData = nan(length(dayBlocks),2);
        for dn = 1:2
            if dn == 1
                dnBlocks = dayBlocks;
            else
                dnBlocks = nightBlocks;
            end

            for ii = 1:length(dnBlocks)
                tmpBlock = dnBlocks{ii};
                [R,P] = corrcoef(smSunful(tmpBlock),nRmACTr(tmpBlock));
                if dn == 1
                    dayData(ii,1) = R(2,1);
                    dayData(ii,2) = P(2,1);
                else
                    nightData(ii,1) = R(2,1);
                    nightData(ii,2) = P(2,1);
                end


            end
        end

        % Interpolate and mean/sd blocks
        for si = 1:2 % daynight blocks
            switch si
                case 1 % day
                    % Find longest block
                    dayBlMax = max(cellfun(@(x) size(x,2), dayBlocks, ...
                        'UniformOutput',true));
                    dayMatlfp = nan(size(dayBlocks,2),dayBlMax);
                    dayMatact = nan(size(dayBlocks,2),dayBlMax);
                    for di = 1:length(dayBlocks)
                        tmpBlock = dayBlocks{di};
                        % LFP
                        lfpDat = smSunful(tmpBlock);
                        % Act
                        actDat = nRmACTr(tmpBlock);
                        oldPoints = 1:length(tmpBlock);
                        newPoints = linspace(1,length(tmpBlock),dayBlMax);

                        newYlfp = interp1(oldPoints,lfpDat,newPoints,'spline');
                        newYact = interp1(oldPoints,actDat,newPoints,'spline');
%                         plot(oldPoints,lfpDat,'o',newPoints,newYlfp,':.');
%                         plot(oldPoints,actDat,'o',newPoints,newYact,':.');
                        dayMatlfp(di,:) = newYlfp;
                        dayMatact(di,:) = newYact;
                    end
                case 2
                    nightBlMax = max(cellfun(@(x) size(x,2), nightBlocks, ...
                        'UniformOutput',true));
                    nightMatlfp = nan(size(nightBlocks,2),nightBlMax);
                    nightMatact = nan(size(nightBlocks,2),nightBlMax);
                    for ni = 1:length(nightBlocks)
                        tmpBlock = nightBlocks{ni};
                        % LFP
                        lfpDat = smSunful(tmpBlock);
                        % Act
                        actDat = nRmACTr(tmpBlock);
                        oldPoints = 1:length(tmpBlock);
                        newPoints = linspace(1,length(tmpBlock),nightBlMax);

                        newYlfp = interp1(oldPoints,lfpDat,newPoints,'spline');
                        newYact = interp1(oldPoints,actDat,newPoints,'spline');
%                         plot(oldPoints,lfpDat,'o',newPoints,newYlfp,':.');
%                         plot(oldPoints,actDat,'o',newPoints,newYact,':.');
                        nightMatlfp(ni,:) = newYlfp;
                        nightMatact(ni,:) = newYact;
                    end
            end
        end

        % Filled significant
        % Unfilled non-signficant
        daySig = dayData(:,2) < 0.5;
        nightSig = nightData(:,2) < 0.5;
        violinDATA = [{nightData(nightSig,1)} {dayData(daySig,1)}];
        catnames_labels = {'Day','Night'};
        coloRS = {cMAP(10,:), cMAP(256,:)};
        violinplot(violinDATA,catnames_labels,'ViolinColor',coloRS,'ViolinAlpha',{0.5 0.5},...
            'ShowMedian',false);
        title('Cross-Cor for night and day epochs')
        xticklabels('Significant correlations')
        ylabel('\textit{r}','Interpreter','latex')

        % Plot 5 ##########################################################
        % Wake up and To bed Event mean plots
        % To Do
        % 3. Add STD
        nexttile(11)

        gbFFT = evData.GoingToBed.FFTBinData;
        gbHz = evData.GoingToBed.Frequency;
        gbFFTt = gbFFT(1:82,:);
        gbHzt = gbHz(1:82,:);
        wuFFT = evData.WakingUp.FFTBinData;
        wuHz = evData.WakingUp.Frequency;
        wuHzt = wuHz(1:82,:);
        wuFFTt = wuFFT(1:82,:);
        % Smooth - individually
        gFFTs = zeros(size(wuFFTt));
        wFFTs = zeros(size(wuFFTt));
        for smo = 1:2
            if smo == 1
                for iii = 1:size(gbFFTt,2)
                    tmpCol = gbFFTt(:,iii);
                    smFFT = smoothdata(tmpCol,'gaussian',6);
                    gFFTs(:,iii) = smFFT;
                end
            else
                for iii = 1:size(wuFFTt,2)
                    tmpCol = wuFFTt(:,iii);
                    smFFT = smoothdata(tmpCol,'gaussian',6);
                    wFFTs(:,iii) = smFFT;
                end
            end
        end


        % Normalize - [upack and repack]
        gwBoth = [gFFTs , wFFTs];
        allNunpk = gwBoth(:);
        allNorm1 = normalize(allNunpk, 'range');
        allNormF = reshape(allNorm1,size(gwBoth));

        gNormF = allNormF(:,1:size(gFFTs,2));
        wNormF = allNormF(:,size(gFFTs,2)+1:end);

        % ttest2(gNormF(:),wNormF(:))

        gMean = mean(gNormF,2);
        wMean = mean(wNormF,2);

        gSTD = std(gNormF,[],2);
        gSTDu = transpose(gMean + (gSTD));
        gSTDd = transpose(gMean - (gSTD));

        wSTD = std(wNormF,[],2);
        wSTDu = transpose(wMean + (wSTD));
        wSTDd = transpose(wMean - (wSTD));

        lp1 = plot(gMean,'LineWidth',2.5);
        lp1.Color = cMAP(246,:);
        hold on

        pch1 = patch([1:length(gMean) fliplr(1:length(gMean))],...
            [gSTDd fliplr(gSTDu)],'k');

        pch1.FaceColor = cMAP(246,:);
        pch1.FaceAlpha = 0.3;
        pch1.EdgeColor = 'none';

        lp2 = plot(wMean,'LineWidth',2.5); % Get into bed
        lp2.Color = cMAP(10,:);

        pch2 = patch([1:length(wMean) fliplr(1:length(wMean))],...
            [wSTDd fliplr(wSTDu)],'k');

        pch2.FaceColor = cMAP(10,:);
        pch2.FaceAlpha = 0.3;
        pch2.EdgeColor = 'none';

        xlt = xline(7,'-','theta');
        xlt.LabelVerticalAlignment = 'middle';
        xlt.LabelHorizontalAlignment = 'center';

        xla = xline(12,'-','alpha');
        xla.LabelVerticalAlignment = 'middle';
        xla.LabelHorizontalAlignment = 'center';

        xlb = xline(29,'-','beta');
        xlb.LabelVerticalAlignment = 'middle';
        xlb.LabelHorizontalAlignment = 'center';

%         xline(8,'theta','LabelVerticalAlignment','middle','LabelOrientation','aligned')
%         xline(12,'alpha','LabelVerticalAlignment','middle','LabelOrientation','aligned')
%         xline(29,'beta','LabelVerticalAlignment','middle','LabelOrientation','aligned')
%         xline(60,'gamma','LabelVerticalAlignment','middle','LabelOrientation','aligned')

        xlim([0 70])
        yticks([0 0.25 0.5 0.75])
        ylim([0 0.75])
        ylabel('Scaled power')
        text(40,0.1,'\textit{p} = 1.18e-07','Interpreter','latex')
        legend('Going to bed: mean','Going to bed: SD',...
               'Wake up: mean','Wake up: SD')
        title('Patient Event Markers')
        xlabel('Frequency (Hz)')

%         nexttile(12)

        xVALS = -0.1:0.001:1;

        [gbTheta , wuTheta] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 't');
        [gbAlpha , wuAlpha] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 'a');
        [gbBeta , wuBeta] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 'b');
        [gbGamma , wuGamma] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 'g');
        figure;
        subplot(4,1,1) % theta
        hold on
        plot(xVALS,gbTheta,'Color',cMAP(246,:),'LineWidth',2)
        plot(xVALS,wuTheta,'Color',cMAP(10,:),'LineWidth',2)
        xlim([0 0.7])
        title('theta')
        subplot(4,1,2)

        hold on
        plot(xVALS,gbAlpha,'Color',cMAP(246,:),'LineWidth',2)
        plot(xVALS,wuAlpha,'Color',cMAP(10,:),'LineWidth',2)
        xlim([0 0.7])
        title('alpha')
        subplot(4,1,3)

        hold on
        plot(xVALS,gbBeta,'Color',cMAP(246,:),'LineWidth',2)
        plot(xVALS,wuBeta,'Color',cMAP(10,:),'LineWidth',2)
        xlim([0 0.7])
        title('beta')
        subplot(4,1,4)

        hold on
        plot(xVALS,gbGamma,'Color',cMAP(246,:),'LineWidth',2)
        plot(xVALS,wuGamma,'Color',cMAP(10,:),'LineWidth',2)
        xlim([0 0.7])
        title('gamma')
        xlabel('Scaled power')
        ylabel('Probability density')






    case 2




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




function [dayBlocks , nightBlocks] = getBlocks(alldatBlocks)

% 0 = night
% 1 = day

% nighDayid = nan(1,50);

nblkC = 1;
nepc = 1;
nallepochs = cell(1,50);

dblkC = 1;
depc = 1;
dallepochs = cell(1,50);
for dn = 1:2 % loop through day and night search

    switch dn
        case 1
            swID = 0;
            for ei = 1:length(alldatBlocks)
                tmpEi = alldatBlocks(ei);
                if isnan(tmpEi) && nepc == 1
                    continue
                elseif isnan(tmpEi) && nepc ~= 1
                    nepc = 1;
                    nblkC = nblkC + 1;
                elseif tmpEi == swID
                    nallepochs{nblkC}(nepc) = ei;
                    nepc = nepc + 1;
                elseif tmpEi ~= swID && nepc == 1
                    continue
                else
                    nepc = 1;
                    nblkC = nblkC + 1;
                end
            end
            nightBlocks  = nallepochs(cellfun(@(x) ~isempty(x), nallepochs,...
                'UniformOutput',true));
        case 2
            swID = 1;
            for ei = 1:length(alldatBlocks)
                tmpEi = alldatBlocks(ei);
                if isnan(tmpEi) && depc == 1
                    continue
                elseif isnan(tmpEi) && depc ~= 1
                    depc = 1;
                    dblkC = dblkC + 1;
                elseif tmpEi == swID
                    dallepochs{dblkC}(depc) = ei;
                    depc = depc + 1;
                elseif tmpEi ~= swID && depc == 1
                    continue
                else
                    depc = 1;
                    dblkC = dblkC + 1;
                end
            end
            dayBlocks  = dallepochs(cellfun(@(x) ~isempty(x), dallepochs,...
                'UniformOutput',true));
    end
end

end




function [dayMean , daySTDud] = getMeanSTD(inMat , offSET)

dayMean = mean(inMat);

daySTD = std(inMat);

daySTDu = dayMean + (daySTD*offSET);
daySTDd = dayMean - (daySTD*offSET);

daySTDud(1,:) = daySTDu;
daySTDud(2,:) = daySTDd;

end



function [freqBand_fitG , freqBand_fitW] = getBandDat(normG , HZg, normW , HZw, bAND)

switch bAND
    case 't'
        low = 4;
        high = 8;
    case 'a'
        low = 9;
        high = 12;
    case 'b'
        low = 13;
        high = 30;
    case 'g'
        low = 31;
        high = 75;
end


xVALS = -0.1:0.001:1;
gbBeta = normG(HZg(:,1) >= low & HZg(:,1) <= high,:);
gbBetaU = gbBeta(:);
wuBeta = normW(HZw(:,1) >= low & HZw(:,1) <= high,:);
wuBetaU = wuBeta(:);

gbB_pdSix = fitdist(gbBetaU,'Kernel','Width',0.05);
freqBand_fitG = pdf(gbB_pdSix,xVALS);
wuB_pdSix = fitdist(wuBetaU,'Kernel','Width',0.05);
freqBand_fitW = pdf(wuB_pdSix,xVALS);

end

