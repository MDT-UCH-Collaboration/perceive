function [] = makeAlexPerceptFigures_v4(figureNUM,subID,hemi,evFlag,dirPRE)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


switch dirPRE
    case 1 % Home Desktop
        dPrefix = 'D:\Dropbox\';
        gPrefix = 'C:\Users\John\';
    case 2 % Work Desktop
        dPrefix = 'C:\Users\Admin\Dropbox\';
        gPrefix = 'C:\Users\Admin\';
    case 3


end


rostLOC = [dPrefix,'Publications_Meta\InProgress\ABaumgartner_Percept2020'];
cd(rostLOC)
rostER = readtable('SubRoster.csv');

mainLOC = [dPrefix,'Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav'];
cd(mainLOC)
addpath([gPrefix,'Documents\GitHub\perceive\AlexB_SleepStudy']);
addpath([gPrefix,'Documents\GitHub\perceive\AlexB_SleepStudy\matplotlib03232022']);
close all

% Create and load summary key



% To do
%. 1 Fix Figure 1A
%. 2 Add x-label to 1E [Frequency Hz]
%. 3 Add xline d t a b g

% 2.	Summary the cross correlation - Summary by night/day block

% 3a.	circe dieum - Summary of time temporal shuffle data - One value per subject
% 3b.   circe dieum - Summary of variance explained - One value per subject

% 3c.   heat plot of all ? [z-scored lfp]
% 6.	Subset of patients with neuroimaging – correlate with location of contact
%
% 1.	Unique plot for 1-2 patients bilateral recording
% 2.    Unique plot for subject with switched Contact


% subID = '3'
% hemi = 'L'

switch figureNUM
    case 1

        % Set up Figure window
        % 4 ROWS | 4 COLUMNS
        % ROW 1 - Cols 1-3 = LFP plot | Col 4 Heat plot
        mainFig = figure;
        set(mainFig,'Position', [1485 311 1128 863]);
        tiledlayout(3,4,"Padding","tight");
        % Load data
        subjectID = subID;
        hemisphere = hemi;
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');

        if evFlag
            [evData] = getPatDat(subjectID , hemisphere , 'Events');
        end

        [apData] = getPatDat(subjectID , hemisphere , 'ActALL');
        [cleanApT] = trim144Apdata(tmData , apData);

        nLFP = tmData.LFP;
        unfurlLFP = nLFP(:);
        mSunful = unfurlLFP - (min(unfurlLFP));
        mSunful(mSunful > 2.2999e+09) = nan;
        mSunful = normalize(mSunful, 'range');
        smSunful = smoothdata(mSunful,'rloess',10,'omitnan');
        firstDayTime8At8p = [find(tmData.hour(:,4) == 8,1,'first')+144 ...
            find(tmData.hour(:,4) == 20,1,'first')+144] ;

        maxVale = max(smSunful);

        % Plot 1 ##########################################################
        nexttile([1 4])

        cMAP = cividis;
        %         [reMAP1] = reMapCmap(smSunful,cMAP,smSunful,0,'median');
        [reMAP] = reMapCmap(tmData,cMAP,smSunful,1,'timeBased');
        lightCM = cMAP(246,:);
        darkCM = cMAP(10,:);
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

        if evFlag
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

        end

        % Add 6AM to 6PM markers
        xl6a = xline(firstDayTime8At8p(1),'--','8AM','LabelVerticalAlignment','bottom',...
            'LabelHorizontalAlignment','center');
        xl6a.Color = cMAP(246,:);
        xl6p = xline(firstDayTime8At8p(2),'--','8PM','LabelVerticalAlignment','bottom',...
            'LabelHorizontalAlignment','center');
        xl6p.Color = cMAP(50,:);

        % Add grey patch code [above and below]
        % First find matching IN and OUT of bed days
        % ASSUMING IN BED prior to MIDNIGHT AND OUT OF BED AFTER MIDNIGHT next DAY
        if evFlag
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

        title(['Patient ',subID])

        % PLOT 2 ##########################################################
        nexttile([1 4])

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
        p1.Color = lightCM;
        p1.Marker = 'o';
        p2 = plot(xAxisAct(swFinMat == 0),zeros(size(xAxisAct(swFinMat == 0))));
        p2.Color = darkCM;
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
        if evFlag
            nexttile(10)
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

                    if length(tmpBlock) == 1
                        continue
                    end


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
            nightNloc = isnan(nightData(:,2));
            nightData = nightData(~nightNloc,:);
            nightBlocks = nightBlocks(~nightNloc);
            dayNloc = isnan(dayData(:,2));
            dayData = dayData(~dayNloc,:);
            dayBlocks = dayBlocks(~dayNloc);

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
            coloRS = {darkCM, lightCM};
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

            if isfield(evData,'GoingToBed')
                gbFFT = evData.GoingToBed.FFTBinData;
                gbHz = evData.GoingToBed.Frequency;
            elseif isfield(evData,'ToBed')
                gbFFT = evData.ToBed.FFTBinData;
                gbHz = evData.ToBed.Frequency;
            end
            gbFFTt = gbFFT(1:82,:);
            gbHzt = gbHz(1:82,:);

            if isfield(evData,'GettingOutOfBed')
                wuFFT = evData.GettingOutOfBed.FFTBinData;
                wuHz = evData.GettingOutOfBed.Frequency;
            elseif isfield(evData,'WakingUp')
                wuFFT = evData.WakingUp.FFTBinData;
                wuHz = evData.WakingUp.Frequency;
            else
                wuFFT = evData.AwakeInMorning.FFTBinData;
                wuHz = evData.AwakeInMorning.Frequency;
            end
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



            if ~matches(subID,'6')

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
                lp1.Color = lightCM;
                hold on

                pch1 = patch([1:length(gMean) fliplr(1:length(gMean))],...
                    [gSTDd fliplr(gSTDu)],'k');

                pch1.FaceColor = lightCM;
                pch1.FaceAlpha = 0.3;
                pch1.EdgeColor = 'none';

                lp2 = plot(wMean,'LineWidth',2.5); % Get into bed
                lp2.Color = darkCM;

                pch2 = patch([1:length(wMean) fliplr(1:length(wMean))],...
                    [wSTDd fliplr(wSTDu)],'k');

                pch2.FaceColor = darkCM;
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

                nexttile(12)
                xVALS = -0.1:0.001:1;

                [gbTheta , wuTheta, gT , wT] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 't');
                [gbAlpha , wuAlpha, gA , wA] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 'a');
                [gbBeta , wuBeta, gB , wB] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 'b');
                [gbGamma , wuGamma, gG , wG] = getBandDat(gNormF , gbHzt, wNormF , wuHzt, 'g');

                % Prep for kruskal wallis / ANOVA
                gALL = [gT ; gA ; gB ; gG];
                wALL = [wT ; wA ; wB ; wG];
                allDATA = [gALL ; wALL];
                bandIDs1 = [repmat({'T'},size(gT)) ; repmat({'A'},size(gA)) ;...
                    repmat({'B'},size(gB)) ; repmat({'G'},size(gG))];
                bandIDs2 = [bandIDs1 ; bandIDs1];
                groupIDs = [repmat({'GO'},size(gALL)) ; repmat({'GET'},size(wALL))];


                %         [p,tbl] = anova2(allDATA,groupIDs,bandIDs2,'off')
                [~,~,stats]  = anovan(allDATA,{bandIDs2,groupIDs},'model',2,'varnames',{'BAND','ToD'},'display','off');
                [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'display','off');

                tbl = array2table(results,"VariableNames", ...
                    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                tbl.("Group A")=gnames(tbl.("Group A"));
                tbl.("Group B")=gnames(tbl.("Group B"));

                % Get max
                maxALL = round(max([gbTheta , wuTheta , gbAlpha , wuAlpha ,...
                    gbBeta, wuBeta,gbGamma , wuGamma]));
                ylinOFF1 = linspace(maxALL,round(maxALL*4),4);
                ylinOFF2 = linspace(0,maxALL*3,4);
                offSET1 = [0 0.5 1 1.5];
                ylinES = ylinOFF1 + offSET1;
                xVALSoff = ylinOFF2 + offSET1;
                hold on
                % PLOT 1
                plot(xVALS,gbTheta + xVALSoff(1),'Color',lightCM,'LineWidth',2)
                plot(xVALS,wuTheta,'Color',darkCM,'LineWidth',2)
                yline(ylinES(1),'-','theta','LabelVerticalAlignment','bottom')
                % Stat
                pVALt = tbl.("P-value")(matches(tbl.("Group A"),'BAND=T,ToD=GO') &...
                    matches(tbl.("Group B"),'BAND=T,ToD=GET'));
                if pVALt == 0
                    ptextT = 'p < 0.0001';
                else
                    pNUM = num2str(round(pVALt,2,'significant'));
                    ptextT = ['p = ',pNUM];
                end
                text(0.3,ylinES(1)-3.5,ptextT)
                % PLOT 2
                plot(xVALS,gbAlpha + xVALSoff(2),'Color',lightCM,'LineWidth',2)
                plot(xVALS,wuAlpha+8.5,'Color',darkCM,'LineWidth',2)
                yline(ylinES(2),'-','alpha','LabelVerticalAlignment','bottom')
                % Stat
                pVALa = tbl.("P-value")(matches(tbl.("Group A"),'BAND=A,ToD=GO') &...
                    matches(tbl.("Group B"),'BAND=A,ToD=GET'));
                if pVALt == 0
                    ptextA = 'p < 0.0001';
                else
                    pNUM = num2str(round(pVALa,2,'significant'));
                    ptextA = ['p = ',pNUM];
                end
                text(0.3,ylinES(2)-3.5,ptextA)
                % PLOT 3
                plot(xVALS,gbBeta + xVALSoff(3),'Color',lightCM,'LineWidth',2)
                plot(xVALS,wuBeta+17,'Color',darkCM,'LineWidth',2)
                yline(ylinES(3),'-','beta','LabelVerticalAlignment','bottom')
                % Stat
                pVALb = tbl.("P-value")(matches(tbl.("Group A"),'BAND=B,ToD=GO') &...
                    matches(tbl.("Group B"),'BAND=B,ToD=GET'));
                if pVALb == 0
                    ptextB = 'p < 0.0001';
                else
                    pNUM = num2str(round(pVALb,2,'significant'));
                    ptextB = ['p = ',pNUM];
                end
                text(0.3,ylinES(3)-3.5,ptextB)
                % PLOT 4
                plot(xVALS,gbGamma + xVALSoff(4),'Color',lightCM,'LineWidth',2)
                plot(xVALS,wuGamma+25.5,'Color',darkCM,'LineWidth',2)
                yline(ylinES(4)+2,'-','gamma','LabelVerticalAlignment','bottom')
                % Stat
                pVALg = tbl.("P-value")(matches(tbl.("Group A"),'BAND=G,ToD=GO') &...
                    matches(tbl.("Group B"),'BAND=G,ToD=GET'));
                if pVALg == 0
                    ptextG = 'p < 0.0001';
                else
                    pNUM = num2str(round(pVALg,2,'significant'));
                    ptextG = ['p = ',pNUM];
                end
                text(0.3,ylinES(4)-3.5,ptextG)

                xlim([0 0.7])
                ylim([0 35.5])
                yticks(linspace(1,7,3))

                xlabel('Scaled power')
                ylabel('Probability density')
                title('Event Frequency bands')

            end

        end

    case 2
        % Make prophet plot
        proph = readtable('SPPD3_L_Forcast.csv');
        % ds x-axis
        % y-hat y axis - dark blue line
        % data y axis - black points
        % y-hat upper/lower - light blue patch


        % SVM - for predicting sleep/wake from LFP

    case 3

        cMAP = cividis;
        lightCM = cMAP(246,:);
        darkCM = cMAP(10,:);
        subjectID = subID;
        hemisphere = hemi;
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');

        if evFlag
            [evData] = getPatDat(subjectID , hemisphere , 'Events');
        end

        [apData] = getPatDat(subjectID , hemisphere , 'ActALL');
        [cleanApT] = trim144Apdata(tmData , apData);

        nLFP = tmData.LFP;
        unfurlLFP = nLFP(:);
        mSunful = unfurlLFP - (min(unfurlLFP));
        mSunful(mSunful > 2.2999e+09) = nan;
        mSunful = normalize(mSunful, 'range');
        smSunful = smoothdata(mSunful,'rloess',10,'omitnan');

        % Method one for creating histograms
        % 1. Use Actigraphy night and day to index into LFP
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

        dayLFPs = smSunful(swFinMat == 1);
        %         allHour = tmData.hour(:);
        %         dayHour = allHour(swFinMat == 1);
        nightLFPs = smSunful(swFinMat == 0);
        %         nightHour = allHour(swFinMat == 0);
        %         dayNightDat = [dayLFPs;nightLFPs];
        %         dayNightID = [repmat({'day'},size(dayLFPs)) ; repmat({'night'},size(nightLFPs))];
        %         dayNightHour = [dayHour ; nightHour];

        %         tbl = table(dayNightDat , dayNightID , dayNightHour,'VariableNames',{'LFP','DN','Hour'});
        %
        %         x = categorical(tbl.DN,{'day','night'});
        %         y = tbl.LFP;
        %         swarmchart(x,y,20,'filled');

        close all
        dayHist = histogram(dayLFPs);

        dayHist.Normalization = "probability";
        dayHist.BinWidth = 0.025;
        dayHist.EdgeColor = 'none';
        dayHist.FaceColor = lightCM;
        dayHist.FaceAlpha = 0.1;

        hold on
        dayStairs = stairs([0 ,dayHist.BinEdges],[0 , dayHist.Values , 0],'k');
        dayStairs.Color = lightCM;

        nightHist = histogram(nightLFPs);
        nightHist.Normalization = "probability";
        nightHist.BinWidth = 0.025;
        nightHist.EdgeColor = 'none';
        nightHist.FaceColor = darkCM;
        nightHist.FaceAlpha = 0.1;
        nightStairs = stairs([0 ,nightHist.BinEdges],[0, nightHist.Values , 0],'k');
        nightStairs.Color = darkCM;

        ylim([0 0.25])
        yticks([0 0.125 0.25])
        ylabel('Prob. density')
        xlim([0 1])
        xticks([0 0.5 1])
        xlabel('Scaled LFP')
        legend({'Awake','','Sleep',''})

        axis square
        %         nightSwarm = swarmchart(nightLFPs);


    case 4
        % Collect data for all patients indicated in CASE 3

        hemiID = zeros(height(rostER),1);
        hemiLab = cell(height(rostER),1);
        subLab = cell(height(rostER),1);
        subID = zeros(height(rostER),1);
        lfpPEak = zeros(height(rostER),1);
        awakeLFP = cell(height(rostER),1);
        asleepLFP = cell(height(rostER),1);
        mdLFPawk = zeros(height(rostER),1);
        for ri = 1:height(rostER)

            tmpSUB = num2str(rostER.subID(ri));
            tmpHEMI = upper(rostER.hemI{ri});
            hemiLab{ri} = tmpHEMI;
            subLab{ri} = tmpSUB;
            [tmData] = getPatDat(tmpSUB , tmpHEMI , 'TimeLine');

            [apData] = getPatDat(tmpSUB , tmpHEMI , 'ActALL');
            [cleanApT] = trim144Apdata(tmData , apData);

            lfpPEak(ri) = tmData.senseFreq;

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            mSunful = unfurlLFP - (min(unfurlLFP));
            mSunful(mSunful > 2.2999e+09) = nan;
            mSunful = normalize(mSunful, 'range');

            mdLFPawk(ri) = median(mSunful,'omitnan');

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

            dayLFPs = mSunful(swFinMat == 1);
            nightLFPs = mSunful(swFinMat == 0);

            awakeLFP{ri} = dayLFPs;
            asleepLFP{ri} = nightLFPs;

            if matches(tmpHEMI,'R')
                hemiID(ri) = 1;
            end
            subID(ri) = rostER.subID(ri);

        end

        % Unpack, normalize, and repack
        allDATA = [];
        allpat = [];
        allAA = {};
        for spI = 1:height(awakeLFP)

            % awake
            allDATA = [allDATA ; awakeLFP{spI}];
            allpat = [allpat ; repmat(spI,size(awakeLFP{spI}))];
            allAA = [allAA ; repmat({'awake'},size(awakeLFP{spI}))];
            % asleep
            allDATA = [allDATA ; asleepLFP{spI}];
            allpat = [allpat ; repmat(spI,size(asleepLFP{spI}))];
            allAA = [allAA ; repmat({'asleep'},size(asleepLFP{spI}))];

        end
        allDATAnorm = normalize(allDATA, 'range');

        % repack
        awakeLFP2 = cell(height(awakeLFP),1);
        asleepLFP2 = cell(height(asleepLFP),1);
        medianAW = zeros(height(awakeLFP),1);
        iqrAW = zeros(height(awakeLFP),1);
        medianAS = zeros(height(awakeLFP),1);
        iqrAS = zeros(height(awakeLFP),1);
        for spI2 = 1:height(awakeLFP)
            awakeTmp = ismember(allpat,spI2) & matches(allAA,'awake');
            asleepTmp = ismember(allpat,spI2) & matches(allAA,'asleep');

            awakeLFP2{spI2} = allDATAnorm(awakeTmp);
            medianAW(spI2) = median(awakeLFP2{spI2},'omitnan');
            iqrAW(spI2) = iqr(awakeLFP2{spI2});
            asleepLFP2{spI2} = allDATAnorm(asleepTmp);
            medianAS(spI2) = median(asleepLFP2{spI2},'omitnan');
            iqrAS(spI2) = iqr(asleepLFP2{spI2});
        end

        % Left first
        %         [~ , hemiSortI] = sort(hemiID);
        %         lhemiSlfpAW =  awakeLFP(hemiSortI);
        %         lhemiSlfpAS =  asleepLFP(hemiSortI);
        %         peakSORTh = lfpPEak(hemiSortI);

        % Sort by max norm power
        [~ , mdPwrI] = sort(mdLFPawk,'descend');

        % Beta peak sort
        %         [~ , peakSortI] = sort(hemiSortI);
        lhPkSlfpAW =  awakeLFP2(mdPwrI);
        lhPkSlfpAS =  asleepLFP2(mdPwrI);
        medAWs = medianAW(mdPwrI);
        medASs = medianAS(mdPwrI);
        hemiSor = hemiLab(mdPwrI);
        subSor = subLab(mdPwrI);
        yAXtk = height(lhPkSlfpAS)*2:-2:1;

        %%%%%% HORIZONTAL SWARM CHART
        %         for spI = 1:height(lhPkSlfpAS)
        %
        %             aSLeepX = lhPkSlfpAS{spI};
        %             aSLeepY = repmat(yAXtk(spI),size(aSLeepX));
        %             aSleepSC = swarmchart(aSLeepX,aSLeepY,10,'filled');
        %             aSleepSC.XJitter = 'none';
        %             aSleepSC.YJitter = 'rand';
        %             aSleepSC.YJitterWidth = 0.4;
        %             aSleepSC.MarkerFaceAlpha = 0.2;
        %             aSleepSC.MarkerEdgeColor = "none";
        %             aSleepSC.MarkerFaceColor = [0 0.447 0.741];
        %             hold on
        %             aAwakeX = lhPkSlfpAW{spI};
        %             aAwakeY = repmat(yAXtk(spI)-1,size(aAwakeX));
        %             aAwakeSC = swarmchart(aAwakeX,aAwakeY,10,'filled');
        %             aAwakeSC.XJitter = 'none';
        %             aAwakeSC.YJitter = 'rand';
        %             aAwakeSC.YJitterWidth = 0.4;
        %             aAwakeSC.MarkerFaceAlpha = 0.2;
        %             aAwakeSC.MarkerEdgeColor = "none";
        %             aAwakeSC.MarkerFaceColor = [0.929 0.694 0.125];
        %
        %         end
        %
        %         plot(medASs,yAXtk,'LineStyle','-','LineWidth',2,'Color',[0 0.447 0.741]);
        %         s1 = scatter(medASs,transpose(yAXtk),100,[0 0.447 0.741],'filled');
        %         s1.MarkerEdgeColor = 'k';
        %         plot(medAWs,yAXtk-1,'LineStyle','-','LineWidth',2,'Color',[0.929 0.694 0.125]);
        %         s2 = scatter(medAWs,transpose(yAXtk)-1,100,[0.929 0.694 0.125],'filled');
        %         s2.MarkerEdgeColor = 'k';
        %         xlim([0 1])
        %         ylim([0 33])
        %
        %         axis square
        %         xticks([0 0.5 1])
        %         xlabel('Scaled LFP')
        %         ylabel('All subjects and hemipsheres')
        %         yticks(2:2:32)
        %         subHemiLabs = cellfun(@(x,y) [x , y], subSor , hemiSor, 'UniformOutput',false);
        %         yticklabels(subHemiLabs)

        %%%%%% DumbBell chart o-----------o asleep to awake median
        %%%%%% Bubble by IQR
        xBCAsleep = medASs;
        yBCAsleep = 1:4:length(mdLFPawk)*4;
        szBCAsleep = iqrAS;

        xBCAwake = medAWs;
        yBCAwake = 1:4:length(mdLFPawk)*4;
        szBCAwake = iqrAW;

        line(transpose([xBCAsleep xBCAwake]), [yBCAsleep ; yBCAsleep],'Color','k')

        hold on

        bc_asleep = scatter(xBCAsleep,yBCAsleep,100,'filled');
        bc_asleep.ColorVariable = [0 0.447 0.741];
        bc_asleep.MarkerFaceColor = [0 0.447 0.741];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0 0.447 0.741];
        bc_asleep.MarkerEdgeAlpha = 1;


        bc_asleep = scatter(xBCAwake,yBCAwake,100,'filled');
        bc_asleep.ColorVariable = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceColor = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0.929 0.694 0.125];
        bc_asleep.MarkerEdgeAlpha = 1;

        xlim([0 1])
        ylim([-1 length(mdLFPawk)*4])

        axis square
        xticks([0 0.5 1])
        xlabel('Scaled LFP')
        ylabel('All subjects and hemipsheres')
        yticks(1:4:length(mdLFPawk)*4)
        subHemiLabs = cellfun(@(x,y) [x , y], subSor , hemiSor, 'UniformOutput',false);
        yticklabels(subHemiLabs)
        legend({'','','','','','','','','',...
            '','','','','','','','Asleep','Awake'})

    case 5

        % Combine - swarmchart in 4 and dot plot in 6
        % Each swarmchart line - should combine data - color

        hemiID = zeros(height(rostER),1);
        hemiLab = cell(height(rostER),1);
        subLab = cell(height(rostER),1);
        subID = zeros(height(rostER),1);
        lfpPEak = zeros(height(rostER),1);
        awakeLFP = cell(height(rostER),1);
        asleepLFP = cell(height(rostER),1);
        mdLFPawk = zeros(height(rostER),1);
        for ri = 1:height(rostER)

            tmpSUB = num2str(rostER.subID(ri));
            tmpHEMI = upper(rostER.hemI{ri});
            hemiLab{ri} = tmpHEMI;
            subLab{ri} = tmpSUB;
            [tmData] = getPatDat(tmpSUB , tmpHEMI , 'TimeLine');

            [apData] = getPatDat(tmpSUB , tmpHEMI , 'ActALL');
            [cleanApT] = trim144Apdata(tmData , apData);

            lfpPEak(ri) = tmData.senseFreq;

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            mSunful = unfurlLFP - (min(unfurlLFP));
            mSunful(mSunful > 2.2999e+09) = nan;
            mSunful = normalize(mSunful, 'range');

            mdLFPawk(ri) = median(mSunful,'omitnan');

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

            dayLFPs = mSunful(swFinMat == 1);
            nightLFPs = mSunful(swFinMat == 0);

            awakeLFP{ri} = dayLFPs;
            asleepLFP{ri} = nightLFPs;

            if matches(tmpHEMI,'R')
                hemiID(ri) = 1;
            end
            subID(ri) = rostER.subID(ri);

        end

        % Unpack, normalize, and repack
        allDATA = [];
        allpat = [];
        allAA = {};
        for spI = 1:height(awakeLFP)

            % awake
            allDATA = [allDATA ; awakeLFP{spI}];
            allpat = [allpat ; repmat(spI,size(awakeLFP{spI}))];
            allAA = [allAA ; repmat({'awake'},size(awakeLFP{spI}))];
            % asleep
            allDATA = [allDATA ; asleepLFP{spI}];
            allpat = [allpat ; repmat(spI,size(asleepLFP{spI}))];
            allAA = [allAA ; repmat({'asleep'},size(asleepLFP{spI}))];

        end
        allDATAnorm = normalize(allDATA, 'range');

        % repack
        awakeLFP2 = cell(height(awakeLFP),1);
        asleepLFP2 = cell(height(asleepLFP),1);
        medianAW = zeros(height(awakeLFP),1);
        iqrAW = zeros(height(awakeLFP),1);
        medianAS = zeros(height(awakeLFP),1);
        iqrAS = zeros(height(awakeLFP),1);
        for spI2 = 1:height(awakeLFP)
            awakeTmp = ismember(allpat,spI2) & matches(allAA,'awake');
            asleepTmp = ismember(allpat,spI2) & matches(allAA,'asleep');

            awakeLFP2{spI2} = allDATAnorm(awakeTmp);
            medianAW(spI2) = median(awakeLFP2{spI2},'omitnan');
            iqrAW(spI2) = iqr(awakeLFP2{spI2});
            asleepLFP2{spI2} = allDATAnorm(asleepTmp);
            medianAS(spI2) = median(asleepLFP2{spI2},'omitnan');
            iqrAS(spI2) = iqr(asleepLFP2{spI2});
        end

        % Left first
        %         [~ , hemiSortI] = sort(hemiID);
        %         lhemiSlfpAW =  awakeLFP(hemiSortI);
        %         lhemiSlfpAS =  asleepLFP(hemiSortI);
        %         peakSORTh = lfpPEak(hemiSortI);

        % Sort by max norm power
        [~ , mdPwrI] = sort(mdLFPawk,'descend');

        % Beta peak sort
        %         [~ , peakSortI] = sort(hemiSortI);
        lhPkSlfpAW =  awakeLFP2(mdPwrI);
        lhPkSlfpAS =  asleepLFP2(mdPwrI);
        medAWs = medianAW(mdPwrI);
        medASs = medianAS(mdPwrI);
        hemiSor = hemiLab(mdPwrI);
        subSor = subLab(mdPwrI);
        yAXtk1 = 1:3:height(lhPkSlfpAS)*3;
        yAXtk3 = 3:3:height(lhPkSlfpAS)*3;

        %%%%%% HORIZONTAL SWARM CHART
        for spI = 1:height(lhPkSlfpAS)

            aSLeepX = lhPkSlfpAS{spI};
            aSLeepY = repmat(yAXtk1(spI),size(aSLeepX));
            aSleepSC = swarmchart(aSLeepX,aSLeepY,10,'filled');
            aSleepSC.XJitter = 'none';
            aSleepSC.YJitter = 'rand';
            aSleepSC.YJitterWidth = 0.4;
            aSleepSC.MarkerFaceAlpha = 0.2;
            aSleepSC.MarkerEdgeColor = "none";
            aSleepSC.MarkerFaceColor = [0 0.447 0.741];
            hold on
            aAwakeX = lhPkSlfpAW{spI};
            aAwakeY = repmat(yAXtk3(spI),size(aAwakeX));
            aAwakeSC = swarmchart(aAwakeX,aAwakeY,10,'filled');
            aAwakeSC.XJitter = 'none';
            aAwakeSC.YJitter = 'rand';
            aAwakeSC.YJitterWidth = 0.4;
            aAwakeSC.MarkerFaceAlpha = 0.2;
            aAwakeSC.MarkerEdgeColor = "none";
            aAwakeSC.MarkerFaceColor = [0.929 0.694 0.125];

        end

        axis square
        xticks([0 0.5 1])
        xlabel('Scaled LFP')
        ylabel('All subjects and hemipsheres')
        yticks(2:3:height(lhPkSlfpAS)*3)
        subHemiLabs = cellfun(@(x,y) [x , y], subSor , hemiSor, 'UniformOutput',false);
        yticklabels(subHemiLabs)

        xBCAsleep = medASs;
        yBCAsleep = yAXtk1;

        xBCAwake = medAWs;
        yBCAwake = yAXtk3;

        line(transpose([xBCAsleep xBCAwake]), [yBCAsleep ; yBCAwake],'Color','k')

        hold on

        bc_asleep = scatter(xBCAsleep,yBCAsleep,100,'filled');
        bc_asleep.ColorVariable = [0 0.447 0.741];
        bc_asleep.MarkerFaceColor = [0 0.447 0.741];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0 0.447 0.741];
        bc_asleep.MarkerEdgeAlpha = 1;


        bc_asleep = scatter(xBCAwake,yBCAwake,100,'filled');
        bc_asleep.ColorVariable = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceColor = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0.929 0.694 0.125];
        bc_asleep.MarkerEdgeAlpha = 1;


    case 6

        % Correlate
        % 2. Difference in median with peak-freq
        % 3. Difference in median with stimulation amplitude
        hemiID = zeros(height(rostER),1);
        hemiLab = cell(height(rostER),1);
        subLab = cell(height(rostER),1);
        subID = zeros(height(rostER),1);
        lfpPEak = zeros(height(rostER),1);
        awakeLFP = cell(height(rostER),1);
        asleepLFP = cell(height(rostER),1);
        mdLFPawk = zeros(height(rostER),1);
        awakeSTM = zeros(height(rostER),1);
        asleepSTM = zeros(height(rostER),1);
        for ri = 1:height(rostER)

            tmpSUB = num2str(rostER.subID(ri));
            tmpHEMI = upper(rostER.hemI{ri});
            hemiLab{ri} = tmpHEMI;
            subLab{ri} = tmpSUB;
            [tmData] = getPatDat(tmpSUB , tmpHEMI , 'TimeLine');

            [apData] = getPatDat(tmpSUB , tmpHEMI , 'ActALL');
            [cleanApT] = trim144Apdata(tmData , apData);

            lfpPEak(ri) = tmData.senseFreq;

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            mSunful = unfurlLFP - (min(unfurlLFP));
            mSunful(mSunful > 2.2999e+09) = nan;
            mSunful = normalize(mSunful, 'range');

            % Stim data
            stimDat = tmData.Stim;
            unfurlSTM = stimDat(:);

            mdLFPawk(ri) = median(mSunful,'omitnan');

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

            dayLFPs = mSunful(swFinMat == 1);
            nightLFPs = mSunful(swFinMat == 0);

            daySTIM = unfurlSTM(swFinMat == 1);
            nightSTIM = unfurlSTM(swFinMat == 0);
            awakeSTM(ri) = mean(daySTIM(daySTIM ~= 0),'omitnan');
            asleepSTM(ri) = mean(nightSTIM(nightSTIM ~= 0),'omitnan');

            awakeLFP{ri} = dayLFPs;
            asleepLFP{ri} = nightLFPs;

            if matches(tmpHEMI,'R')
                hemiID(ri) = 1;
            end
            subID(ri) = rostER.subID(ri);

        end

        % Unpack, normalize, and repack
        allDATA = [];
        allpat = [];
        allAA = {};
        for spI = 1:height(awakeLFP)

            % awake
            allDATA = [allDATA ; awakeLFP{spI}];
            allpat = [allpat ; repmat(spI,size(awakeLFP{spI}))];
            allAA = [allAA ; repmat({'awake'},size(awakeLFP{spI}))];
            % asleep
            allDATA = [allDATA ; asleepLFP{spI}];
            allpat = [allpat ; repmat(spI,size(asleepLFP{spI}))];
            allAA = [allAA ; repmat({'asleep'},size(asleepLFP{spI}))];

        end
        allDATAnorm = normalize(allDATA, 'range');

        % repack
        awakeLFP2 = cell(height(awakeLFP),1);
        asleepLFP2 = cell(height(asleepLFP),1);
        medianAW = zeros(height(awakeLFP),1);
        iqrAW = zeros(height(awakeLFP),1);
        medianAS = zeros(height(awakeLFP),1);
        iqrAS = zeros(height(awakeLFP),1);
        for spI2 = 1:height(awakeLFP)
            awakeTmp = ismember(allpat,spI2) & matches(allAA,'awake');
            asleepTmp = ismember(allpat,spI2) & matches(allAA,'asleep');

            awakeLFP2{spI2} = allDATAnorm(awakeTmp);
            medianAW(spI2) = median(awakeLFP2{spI2},'omitnan');
            iqrAW(spI2) = iqr(awakeLFP2{spI2});
            asleepLFP2{spI2} = allDATAnorm(asleepTmp);
            medianAS(spI2) = median(asleepLFP2{spI2},'omitnan');
            iqrAS(spI2) = iqr(asleepLFP2{spI2});
        end

        medOffset = medianAW - medianAS;

        % Get lag
        % Cross Correlation
        % Loop through all patients

        allxcor = nan(height(rostER),1);
        for ri = 1:height(rostER)

            tmpSUB = num2str(rostER.subID(ri));
            tmpHEMI = upper(rostER.hemI{ri});
            hemiLab{ri} = tmpHEMI;
            subLab{ri} = tmpSUB;
            [tmData] = getPatDat(tmpSUB , tmpHEMI , 'TimeLine');
            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            mSunful = unfurlLFP - (min(unfurlLFP));
            mSunful(mSunful > 2.2999e+09) = nan;
            mSunful = normalize(mSunful, 'range');

            [apData] = getPatDat(tmpSUB , tmpHEMI , 'ActALL');
            [cleanApT] = trim144Apdata(tmData , apData);
            nActraw = cleanApT.Activity;
            unfurlACTr = nActraw(:);
            sMnRmACTr = smoothdata(unfurlACTr,'gaussian',40,'omitnan');
            nRmACTr = normalize(sMnRmACTr, 'range');

            nonNanLocs = ~isnan(nRmACTr);

            [crossCorvals,lags] = xcorr(mSunful(nonNanLocs),nRmACTr(nonNanLocs));
            crossCorvalsN = crossCorvals/max(crossCorvals);
            [~,ccnLoc] = max(crossCorvalsN);
            maxClagLoc = lags(ccnLoc);
            if isempty(maxClagLoc)
                allxcor(ri) = nan;
            else
                allxcor(ri) = maxClagLoc;
            end
        end

        %         plot(lags,crossCorvalsN,'Color','k','LineWidth',1.5)

        % TOTAL
        % Scale or normalize first
        allxcor2 = abs(allxcor);
        % get rid of outlier
        remOUTallxco = allxcor2(allxcor < 500);
        remmedOffset = medOffset(allxcor < 500);

        % Get coefficients of a line fit through the data.
        coefficients = polyfit(remmedOffset, remOUTallxco, 1);
        % Create a new x axis with exactly 1000 points (or whatever you want).
        xFit = linspace(min(remmedOffset), max(remmedOffset), 1000);
        % Get the estimated yFit value for each of those 1000 new x locations.
        yFit = polyval(coefficients , xFit);

        dp = plot(remmedOffset,remOUTallxco,'o');
        dp.MarkerFaceColor = 'k';
        dp.MarkerEdgeColor = 'k';
        dp.MarkerSize = 10;
        hold on
        plot(xFit, yFit, 'k--', 'LineWidth', 2);
        ylim([-4 140])
        yticks([0 70 140])
        ylabel('Sample offset between actigraphy and LFP')
        xlim([0 0.35])
        xticks([0 0.175 0.35])
        xlabel('Scaled LFP delta between Asleep and Awake')
        [r,p] = corr(remmedOffset,remOUTallxco);
        text(0.28, 130, ['p = ',num2str(round(p,3))])
        text(0.28, 124, ['r = ',num2str(round(r,3))])
        title('Correlation between median difference in LFP and Lag')

        axis square

        coefficients2 = polyfit(medOffset, lfpPEak, 1);
        % Create a new x axis with exactly 1000 points (or whatever you want).
        xFit2 = linspace(min(medOffset), max(medOffset), 1000);
        % Get the estimated yFit value for each of those 1000 new x locations.
        yFit2 = polyval(coefficients2 , xFit2);

        figure;
        dp2 = plot(medOffset,lfpPEak,'o');
        dp2.MarkerFaceColor = 'k';
        dp2.MarkerEdgeColor = 'k';
        dp2.MarkerSize = 10;
        hold on
        plot(xFit2, yFit2, 'k--', 'LineWidth', 2);
        ylim([0 30])
        yticks([0 15 30])
        ylabel('Peak frequency (Hz)')
        xlim([0 0.35])
        xticks([0 0.175 0.35])
        xlabel('Scaled LFP delta between Asleep and Awake')
        [r2,p2] = corr(medOffset,lfpPEak);
        text(0.28, 28, ['p = ',num2str(round(p2,3))])
        text(0.28, 27, ['r = ',num2str(round(r2,3))])
        title('Correlation between median difference in LFP and Peak Frequency')

        axis square

        allSTM = mean([awakeSTM , asleepSTM],2);
        % Correlate with STIM
        % 1. Correlate Diff with STIM
        coefficients3 = polyfit(medOffset, allSTM, 1);
        % Create a new x axis with exactly 1000 points (or whatever you want).
        xFit3 = linspace(min(medOffset), max(medOffset), 1000);
        % Get the estimated yFit value for each of those 1000 new x locations.
        yFit3 = polyval(coefficients3 , xFit3);


        figure;
        dp2 = plot(medOffset,allSTM,'o');
        dp2.MarkerFaceColor = 'k';
        dp2.MarkerEdgeColor = 'k';
        dp2.MarkerSize = 10;
        hold on
        plot(xFit3, yFit3, 'k--', 'LineWidth', 2);
        ylim([0 7])
        yticks([0 3 6])
        ylabel('Average stimulation (mA)')
        xlim([0 0.35])
        xticks([0 0.175 0.35])
        xlabel('Scaled LFP delta between Asleep and Awake')
        [r3,p3] = corr(allSTM,medOffset);
        text(0.28, 6, ['p = ',num2str(round(p3,3))])
        text(0.28, 5.6, ['r = ',num2str(round(r3,3))])
        title('Correlation between median difference in LFP and stimulation')

        axis square




    case 7 % Cross correlation
        % INDIVIDUAL SESSION? for Cross Cor Analysis
        hemiID = zeros(height(rostER),1);
        hemiLab = cell(height(rostER),1);
        subLab = cell(height(rostER),1);
        subID = zeros(height(rostER),1);
        lfpPEak = zeros(height(rostER),1);
        awakeLFP = cell(height(rostER),1);
        asleepLFP = cell(height(rostER),1);
        mdLFPawk = zeros(height(rostER),1);
        for ri = 1:height(rostER)

            tmpSUB = num2str(rostER.subID(ri));
            tmpHEMI = upper(rostER.hemI{ri});
            hemiLab{ri} = tmpHEMI;
            subLab{ri} = tmpSUB;
            [tmData] = getPatDat(tmpSUB , tmpHEMI , 'TimeLine');

            [apData] = getPatDat(tmpSUB , tmpHEMI , 'ActALL');
            [cleanApT] = trim144Apdata(tmData , apData);

            lfpPEak(ri) = tmData.senseFreq;

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            mSunful = unfurlLFP - (min(unfurlLFP));
            mSunful(mSunful > 2.2999e+09) = nan;
            mSunful = normalize(mSunful, 'range');

            mdLFPawk(ri) = median(mSunful,'omitnan');

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

            dayLFPs = mSunful(swFinMat == 1);
            nightLFPs = mSunful(swFinMat == 0);

            awakeLFP{ri} = dayLFPs;
            asleepLFP{ri} = nightLFPs;

            if matches(tmpHEMI,'R')
                hemiID(ri) = 1;
            end
            subID(ri) = rostER.subID(ri);

        end

        test = 1;






    case 8 % Histogram overlap

        hemiLab = cell(height(rostER),1);
        subLab = cell(height(rostER),1);
        perOVER = nan(height(rostER),1);
        % Loop through subjects / hemispheres
        for ri = 1:height(rostER)

            tmpSUB = num2str(rostER.subID(ri));
            tmpHEMI = upper(rostER.hemI{ri});
            hemiLab{ri} = tmpHEMI;
            subLab{ri} = tmpSUB;
            [tmData] = getPatDat(tmpSUB , tmpHEMI , 'TimeLine');
            [apData] = getPatDat(tmpSUB , tmpHEMI , 'ActALL');
            [cleanApT] = trim144Apdata(tmData , apData);

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            mSunful = unfurlLFP - (min(unfurlLFP));
            mSunful(mSunful > 2.2999e+09) = nan;
            mSunful = normalize(mSunful, 'range');
            smSunful = smoothdata(mSunful,'rloess',10,'omitnan');

            % Method one for creating histograms
            % 1. Use Actigraphy night and day to index into LFP
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

            dayLFPs = smSunful(swFinMat == 1);
            nightLFPs = smSunful(swFinMat == 0);

            bothHistograms = [nightLFPs; dayLFPs];
            minCounts = min(bothHistograms, [], 1);
            maxCounts = max(bothHistograms, [], 1);
            perOVER(ri) = minCounts ./ maxCounts;

        end

        % PLOT STUFF
        absPERO = abs(perOVER);
        [perOVsort , sIND] = sort(absPERO,'ascend');
        hemiLsort = hemiLab(sIND);
        subLsort = subLab(sIND);

        scatter(perOVsort,1:length(perOVsort),50,"black","filled");
        hold on
        scatter(perOVsort(12), 12, 50, "red","filled")
        xlim([0 0.5])
        xticks([0 0.25 0.5])
        xlabel("Fraction overlap between Night/Day LFP")
        ylabel("Subject | Hemisphere ID")
        ylim([0 17])
        yticks(1:length(perOVsort))
        hemiSubL = cellfun(@(x,y) [x , y], subLsort, hemiLsort,'UniformOutput',false);
        yticklabels(hemiSubL);

        axis square


    case 9 % CARPE DIEM

        my_favourite_colour     = [0.8500 0.3250 0.0980];
        subjectID = subID;
        hemisphere = hemi;
        [tmData] = getPatDat(subjectID , hemisphere , 'TimeLine');

        time_stamps = tmData.actTime(:);

        nLFP = tmData.LFP;
        unfurlLFP = nLFP(:);
        mSunful = unfurlLFP - (min(unfurlLFP));
        mSunful(mSunful > 2.2999e+09) = nan;

        values = mSunful;

        figure
        set(gcf,'Units','Normalized','Position',[.2 .4 .6 .25])
        plot_zscored_timeseries(time_stamps, values, my_favourite_colour)

        figure
        set(gcf,'Units','Normalized','Position',[.3 .3 .4 .3])

        time_res        = 1; % time resolution (in hours) for periodogram
        max_period      = 72; % Maximum period of 1 week = 168 hours
        do_normalise    = true; % Whether to normalise the periodogram

        % Calculate periodogram
        [psd_estimate, time_periods] = circadian_periodogram(time_stamps, values, time_res, max_period);

        % Plot the periodogram
        plot_periodogram(psd_estimate, time_periods, do_normalise, my_favourite_colour)

        time_res        = 1;
        n_shuffles      = 200;
        shuffle_type    = 'circshift';

        % Get the proportion of variance explained by time of day
        var_explained = variance_explained_by_timeofday(time_stamps, values, time_res);

        % Get the variance explained for shuffled data n_shuffles times to see whether var_explained is significant
        [~, var_explained_p] = get_shuffled_var_explained(time_stamps, values, time_res, n_shuffles, shuffle_type);

        % Plot a fit based on time of day to the data across days
        figure
        plot_timeofday_fit(time_stamps, values, time_res, my_favourite_colour)

        % Add the variance explained by time of day & p-val to the figure title
        title(['Var explained by TOD: ' num2str(var_explained) ', p =' num2str(var_explained_p)])

        figure
        plot_timeofday_fit(time_stamps, values, time_res, my_favourite_colour,'polar')

        time_res    = 1;
        stat        = 'mean';

        figure
        circadian_rose(time_stamps, values, time_res, stat, my_favourite_colour)

        time_res            = 1; % Temporal resolution of time bins (= heatmap pixels)
        percentile_cutoff   = 2; % For the colour scale of the heatmap, ignore the top and bottom x% of data

        [circadian_matrix, ~] = make_circadian_matrix(time_stamps, values, time_res);

        figure
        plot_circadian_matrix(circadian_matrix, percentile_cutoff, my_favourite_colour);

        close all

        circadian_summary_figure(time_stamps, values,time_res)

        %         % Get shuffled data points
        %         shuffled_data_points    = within_day_shuffle(time_stamps, values, 'circshift');
        %
        %         [circadian_matrix, time_edges] = make_circadian_matrix(time_stamps, shuffled_data_points, time_res);
        %
        %
        %         figure
        %         plot_circadian_matrix(circadian_matrix, percentile_cutoff, my_favourite_colour);



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

            %             tmLoc = find(timMatf == tmpDat);
            [~,tmLoc] = min(abs(timMatf - tmpDat));

            cmapLoc = round((tmLoc/144)*256);

            tCmp = cMap(cmapLoc,:);
            reMAP(di,:) = tCmp;

        end

end

end



function [outIND , dayOUTall] = findBEDinds(sTATE , evDAT , tlDAT)

% check for repeats ***************************************

if matches(sTATE,'inbed')
    if isfield(evDAT,'GoingToBed')
        evDAY = evDAT.GoingToBed.Day;
        evHOUR = evDAT.GoingToBed.Hour;
        evMIN = evDAT.GoingToBed.Minute;
    else
        evDAY = evDAT.ToBed.Day;
        evHOUR = evDAT.ToBed.Hour;
        evMIN = evDAT.ToBed.Minute;
    end
else
    if isfield(evDAT,'GettingOutOfBed')
        evDAY = evDAT.GettingOutOfBed.Day;
        evHOUR = evDAT.GettingOutOfBed.Hour;
        evMIN = evDAT.GettingOutOfBed.Minute;
    elseif isfield(evDAT,'WakingUp')
        evDAY = evDAT.WakingUp.Day;
        evHOUR = evDAT.WakingUp.Hour;
        evMIN = evDAT.WakingUp.Minute;
    else
        evDAY = evDAT.AwakeInMorning.Day;
        evHOUR = evDAT.AwakeInMorning.Hour;
        evMIN = evDAT.AwakeInMorning.Minute;
    end
end

evDAYc = unique(evDAY);
if length(evDAY) ~= length(evDAYc)

    evDAYn = zeros(length(evDAYc),1);
    evHOURn = zeros(length(evDAYc),1);
    evMINn = zeros(length(evDAYc),1);
    for ci = 1:length(evDAYc)
        tmpCD = evDAYc(ci);
        tmpDI = find(tmpCD == evDAY,1,'first');
        evDAYn(ci) = tmpCD;
        evHOURn(ci) = evHOUR(tmpDI);
        evMINn(ci) = evMIN(tmpDI);
    end
    evDAY = evDAYn;
    evHOUR = evHOURn;
    evMIN = evMINn;
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
    % If day doesn't match then skip
    if sum(tldayFind) == 0
        continue
    else
        tlhourFind = ismember(tlHOUR,tmpHOUR) & tldayFind;
        dayOUTall(ei) = tmpDAY;
        tlminFind = find(ismember(tlMIN,tmpMIN) & tlhourFind);
        if isempty(tlminFind)
            hourINDS = find(tlhourFind);
            hourTRIM = tlMIN(tlhourFind);
            [~ , minLOC] = min(abs(hourTRIM - tmpMIN));
            tlminFind = hourINDS(minLOC);
        end
        if isempty(tlminFind)
            continue
        else
            outIND(ei) = tlminFind;
        end
    end


end
% clean up outIND and dayOUTall
outIND = outIND(outIND ~= 0);
dayOUTall = dayOUTall(dayOUTall ~= 0);


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



function [freqBand_fitG , freqBand_fitW , gRAW , wRAW] = getBandDat(normG , HZg, normW , HZw, bAND)

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

gRAW = gbBetaU;
wRAW = wuBetaU;

end

