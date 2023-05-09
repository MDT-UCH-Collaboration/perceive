function [] = alexPerceptR1_npj(dirPRE , plotID)




switch dirPRE
    case 1 % Home Desktop
        dPrefix = 'E:\Dropbox\';
        gPrefix = 'C:\Users\Admin\';
    case 2 % Work Desktop
        dPrefix = 'D:\Dropbox\';
        gPrefix = 'C:\Users\Admin\';
    case 3
        dPrefix = 'C:\Users\johna\Dropbox\';
        gPrefix = 'C:\Users\johna\';

end

%%% FOR NEW FIGURE 2 to depict alpha vs beta
%%% COMBINE 2 and 4

% makeAlexPerceptFigures_v6(3,'3','L',1,3)
rostLOC = [dPrefix,'Publications_Meta\InProgress\ABaumgartner_Percept2020'];
cd(rostLOC)
rostER = readtable('SubRoster.csv');
demoDATA = readtable('demoDATA_12222022.xlsx');

mainLOC = [dPrefix,'Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav'];
dissLOC = [dPrefix,'\Publications_Meta\InProgress\ABaumgartner_Percept2020\DisucssionFigure\outDATA'];
cd(mainLOC)
addpath([gPrefix,'Documents\GitHub\perceive\AlexB_SleepStudy']);
addpath([gPrefix,'Documents\GitHub\perceive\AlexB_SleepStudy\matplotlib03232022']);

revLOC = [dPrefix,'Publications_Meta\InProgress\ABaumgartner_Percept2020\revisionData'];

close all

switch plotID
    case 1


        cd(revLOC)
        load("NewFigureSwarm.mat","lhPkSlfpAS","lhPkSlfpAW","yAXtk1",...
            "yAXtk3","subSor","hemiSor","medASs","medAWs")

        %%%%%% HORIZONTAL SWARM CHART
        for spI = 1:height(lhPkSlfpAS)

            if rostER.eventUSE(spI)
                subplot(1,3,1)

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

            else
                subplot(1,3,2)
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

        end

        subHemiLabs = cellfun(@(x,y) [x , y], subSor , hemiSor, 'UniformOutput',false);
        artLabs = subHemiLabs(~rostER.eventUSE);
        nartLabs = subHemiLabs(logical(rostER.eventUSE));

        subplot(1,3,1)
        axis square
        xticks([0 0.5 1])
        xlabel('Scaled LFP')
        ylabel('All subjects and hemipsheres')
        yticks(2:3:height(lhPkSlfpAS)*3)
        ylim([0 height(lhPkSlfpAS)*3])
        title('Non-artifact events')
        yticks(yAXtk1(logical(rostER.eventUSE)))
        yticklabels(nartLabs)

        subplot(1,3,2)
        axis square
        xticks([0 0.5 1])
        xlabel('Scaled LFP')
        ylabel('All subjects and hemipsheres')
        yticks(2:3:height(lhPkSlfpAS)*3)
        ylim([0 height(lhPkSlfpAS)*3])
        title('Artifact events')
        yticks(yAXtk1(~rostER.eventUSE))
        yticklabels(artLabs)

        xBCAsleepNA = medASs(logical(rostER.eventUSE));
        yBCAsleepNA = yAXtk1(logical(rostER.eventUSE));

        xBCAsleepA = medASs(logical(~rostER.eventUSE));
        yBCAsleepA = yAXtk1(logical(~rostER.eventUSE));

        xBCAwakeNA = medAWs(logical(rostER.eventUSE));
        yBCAwakeNA = yAXtk3(logical(rostER.eventUSE));

        xBCAwakeA = medAWs(logical(~rostER.eventUSE));
        yBCAwakeA = yAXtk3(logical(~rostER.eventUSE));

        subplot(1,3,1)
        line(transpose([xBCAsleepNA xBCAwakeNA]), [yBCAsleepNA ; yBCAwakeNA],'Color','k')

        subplot(1,3,2)
        line(transpose([xBCAsleepA xBCAwakeA]), [yBCAsleepA ; yBCAwakeA],'Color','k')

        hold on

        subplot(1,3,1)
        bc_asleep = scatter(xBCAsleepNA,yBCAsleepNA,70,'filled');
        bc_asleep.ColorVariable = [0 0.447 0.741];
        bc_asleep.MarkerFaceColor = [0 0.447 0.741];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0 0.447 0.741];
        bc_asleep.MarkerEdgeAlpha = 1;
        bc_asleep = scatter(xBCAwakeNA,yBCAwakeNA,70,'filled');
        bc_asleep.ColorVariable = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceColor = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0.929 0.694 0.125];
        bc_asleep.MarkerEdgeAlpha = 1;

        subplot(1,3,2)
        bc_asleep = scatter(xBCAsleepA,yBCAsleepA,70,'filled');
        bc_asleep.ColorVariable = [0 0.447 0.741];
        bc_asleep.MarkerFaceColor = [0 0.447 0.741];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0 0.447 0.741];
        bc_asleep.MarkerEdgeAlpha = 1;
        bc_asleep = scatter(xBCAwakeA,yBCAwakeA,70,'filled');
        bc_asleep.ColorVariable = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceColor = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0.929 0.694 0.125];
        bc_asleep.MarkerEdgeAlpha = 1;


        subplot(1,3,3)
        % Add mean difference bar chart comparison
        percentNA = (xBCAwakeNA - xBCAsleepNA)./xBCAwakeNA;
        percentA = (xBCAwakeA - xBCAsleepA)./xBCAwakeA;
        boxLABels = [ones(length(percentNA),1);...
            ones(length(percentA),1)+1];
        allDATA = [percentNA ; percentA];
        dataTable = table(boxLABels, allDATA, 'VariableNames',{'Group','Data'});

        boxchart(dataTable.Group,dataTable.Data);
        xticks([1 2])
        xticklabels({'Non-Artifact','Artifact'})
        ylabel('Percent difference between awake and asleep LFP')
        yticks([0 0.35 0.7])
        axis square

        [a,b,c] = ranksum(percentNA , percentA)


    case 2


        cd(revLOC)
        load("NewFigureSwarm.mat","lhPkSlfpAS","lhPkSlfpAW","yAXtk1",...
            "yAXtk3","subSor","hemiSor","medASs","medAWs","lfpPEak")

        rostER.lfpPEAK = lfpPEak;

        newTAB = table(subSor, hemiSor);
        for si = 1:height(newTAB)

            tmpSUB = newTAB.subSor{si};
            tmpHEM = newTAB.hemiSor{si};

            rowINDEX = ismember(rostER.subID,str2double(tmpSUB)) &...
                ismember(rostER.hemI,lower(tmpHEM));

            newTAB.LFPpeak(si) = rostER.lfpPEAK(rowINDEX);

        end

        %%%%%% HORIZONTAL SWARM CHART
        for spI = 1:height(lhPkSlfpAS)

            if newTAB.LFPpeak(spI) < 13
                subplot(1,3,1)

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

            else
                subplot(1,3,2)
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

        end

        subHemiLabs = cellfun(@(x,y) [x , y], subSor , hemiSor, 'UniformOutput',false);
        alLabs = subHemiLabs(newTAB.LFPpeak < 13);
        beLabs = subHemiLabs(newTAB.LFPpeak > 13);

        subplot(1,3,1)
        axis square
        xticks([0 0.5 1])
        xlabel('Scaled LFP')
        ylabel('All subjects and hemipsheres')
        yticks(2:3:height(lhPkSlfpAS)*3)
        ylim([0 height(lhPkSlfpAS)*3])
        title('Alpha peaks')
        yticks(yAXtk1(newTAB.LFPpeak < 13))
        yticklabels(alLabs)

        subplot(1,3,2)
        axis square
        xticks([0 0.5 1])
        xlabel('Scaled LFP')
        ylabel('All subjects and hemipsheres')
        yticks(2:3:height(lhPkSlfpAS)*3)
        ylim([0 height(lhPkSlfpAS)*3])
        title('Beta peaks')
        yticks(yAXtk1(newTAB.LFPpeak > 13))
        yticklabels(beLabs)



        xBCAsleepA = medASs(newTAB.LFPpeak < 13);
        yBCAsleepA = yAXtk1(newTAB.LFPpeak < 13);

        xBCAsleepB = medASs(newTAB.LFPpeak > 13);
        yBCAsleepB = yAXtk1(newTAB.LFPpeak > 13);

        xBCAwakeA = medAWs(newTAB.LFPpeak < 13);
        yBCAwakeA = yAXtk3(newTAB.LFPpeak < 13);

        xBCAwakeB = medAWs(newTAB.LFPpeak > 13);
        yBCAwakeB = yAXtk3(newTAB.LFPpeak > 13);

        subplot(1,3,1)
        line(transpose([xBCAsleepA xBCAwakeA]), [yBCAsleepA ; yBCAwakeA],'Color','k')

        subplot(1,3,2)
        line(transpose([xBCAsleepB xBCAwakeB]), [yBCAsleepB ; yBCAwakeB],'Color','k')

        hold on

        subplot(1,3,1)
        bc_asleep = scatter(xBCAsleepA,yBCAsleepA,70,'filled');
        bc_asleep.ColorVariable = [0 0.447 0.741];
        bc_asleep.MarkerFaceColor = [0 0.447 0.741];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0 0.447 0.741];
        bc_asleep.MarkerEdgeAlpha = 1;
        bc_asleep = scatter(xBCAwakeA,yBCAwakeA,70,'filled');
        bc_asleep.ColorVariable = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceColor = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0.929 0.694 0.125];
        bc_asleep.MarkerEdgeAlpha = 1;

        subplot(1,3,2)
        bc_asleep = scatter(xBCAsleepB,yBCAsleepB,70,'filled');
        bc_asleep.ColorVariable = [0 0.447 0.741];
        bc_asleep.MarkerFaceColor = [0 0.447 0.741];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0 0.447 0.741];
        bc_asleep.MarkerEdgeAlpha = 1;
        bc_asleep = scatter(xBCAwakeB,yBCAwakeB,70,'filled');
        bc_asleep.ColorVariable = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceColor = [0.929 0.694 0.125];
        bc_asleep.MarkerFaceAlpha = 1;
        bc_asleep.MarkerEdgeColor = [0.929 0.694 0.125];
        bc_asleep.MarkerEdgeAlpha = 1;

        subplot(1,3,3)
        % Add mean difference bar chart comparison
        percentAl = (xBCAwakeA - xBCAsleepA)./xBCAwakeA;
        percentBe = (xBCAwakeB - xBCAsleepB)./xBCAwakeB;
        boxLABels = [ones(length(percentAl),1);...
            ones(length(percentBe),1)+1];
        allDATA = [percentAl ; percentBe];
        dataTable = table(boxLABels, allDATA, 'VariableNames',{'Group','Data'});

        boxchart(dataTable.Group,dataTable.Data);
        xticks([1 2])
        xticklabels({'Non-Artifact','Artifact'})
        ylabel('Percent difference between awake and asleep LFP')
        yticks([0 0.35 0.7])
        axis square

        [a,b,c] = ranksum(percentAl , percentBe)



    case 3


        all2bd = nan(height(rostER),1);
        allwo = nan(height(rostER),1);

        for ri = 1:height(rostER)
            close all
            tmpSUB = num2str(rostER.subID(ri));
            tmpHEMI = upper(rostER.hemI{ri});
            [evData] = getPatDat(tmpSUB , tmpHEMI , 'Events');

            if ~isstruct(evData)
                continue
            end

            if isfield(evData,'GoingToBed')
                num2bed = length(evData.GoingToBed.DateTime);
            elseif isfield(evData,'ToBed')
                num2bed = length(evData.ToBed.DateTime);
            end
            all2bd(ri) = num2bed;


            if isfield(evData,'GettingOutOfBed')
                numoutB = length(evData.GettingOutOfBed.DateTime);
            elseif isfield(evData,'WakingUp')
                numoutB = length(evData.WakingUp.DateTime);
            else
                numoutB = length(evData.AwakeInMorning.DateTime);
            end
            allwo(ri) = numoutB;
        end


        averageINbed = mean(all2bd,'omitnan');
        averageOUTbed = mean(allwo ,'omitnan');
        sumINbed = sum(all2bd ,'omitnan');
        sumOUTbed = sum(allwo ,'omitnan');
        stdINbed = std(all2bd,'omitnan');
        stdOUTbed = std(allwo,'omitnan');


    case 4


        mainFig = figure;
        set(mainFig,'Position', [1485 311 1128 863]);
        tiledlayout(2,4,"Padding","tight");
        % Load data
        % GOOD SUBJECt
        % subjectID = '3';
        % hemisphere = 'L';


        % ARTIFACT 1
        % subjectID = '9';
        % hemisphere = 'L';

        % ARTIFACT 2
        % subjectID = '3';
        % hemisphere = 'L';

        % Non-ARTIFACT 1
        % subjectID = '8';
        % hemisphere = 'L';

        % Non-ARTIFACT 2
        % subjectID = '4';
        % hemisphere = 'R';

        subjectID = {'9','3','8','4'};
        hemisphere = {'L','L','L','R'};
        titlEE = {'Artifact','Artifact','Non-Artifact','Non-Artifact'};
        for ii = 1:4


            [tmData] = getPatDat(subjectID{ii} , hemisphere{ii} , 'TimeLine');

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            raw_LFP = unfurlLFP - (min(unfurlLFP));
            clean_LFP = raw_LFP;
            clean_LFP(clean_LFP > 2.2999e+09) = nan;
            % clean_LFP(clean_LFP < 50) = nan;
            mSunful = normalize(clean_LFP, 'range');
            smSunful = smoothdata(mSunful,'rloess',10,'omitnan');

            % Plot 1 ##########################################################
            % Raw
            % nexttile([1 4])
            %
            % s1 = scatter(1:length(raw_LFP),raw_LFP,[],'k','filled');
            % s1.MarkerFaceAlpha = 0.5;
            %
            % maxVale = max(raw_LFP);
            % ylim([0 round(maxVale + 0.1,1)])
            % yticks([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
            % yticklabels([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
            % ylabel('Raw LFP power')
            % dayStarts = round(linspace(1,length(raw_LFP)-144,length(raw_LFP)/144));
            % xticks(dayStarts);
            % xticklabels(1:length(raw_LFP)/144)
            % xlim([1, length(raw_LFP)])
            % xlabel('Days of recording')
            % set(gca,'TickLength',[0 .001])

            % Plot 2 ##########################################################
            % Outliers removed
            % nexttile([1 4])

            % s2 = scatter(1:length(clean_LFP),clean_LFP,[],'r','filled');
            % s2.MarkerFaceAlpha = 0.5;
            %
            % numRemPts = sum(isnan(clean_LFP));
            % subtitle(['Number of removed points = ',num2str(numRemPts)],'FontSize',10,'Color','r');
            % ax = gca;
            % ax.TitleHorizontalAlignment = 'left';
            %
            % maxVale = max(clean_LFP);
            % ylim([0 round(maxVale + 0.1,1)])
            % yticks([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
            % yticklabels([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
            % ylabel('Raw LFP power')
            % dayStarts = round(linspace(1,length(clean_LFP)-144,length(clean_LFP)/144));
            % xticks(dayStarts);
            % xticklabels(1:length(clean_LFP)/144)
            % xlim([1, length(clean_LFP)])
            % xlabel('Days of recording')
            % set(gca,'TickLength',[0 .001])

            % Plot 3 ##########################################################
            % Clean and prettified
            nexttile([1 2])
            cMAP = cividis;
            [reMAP] = reMapCmap(tmData,cMAP,smSunful,1,'timeBased');
            scatter(1:length(smSunful),smSunful,[],reMAP,'filled')

            maxVale = max(smSunful);
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
            title([subjectID{ii},'  ',hemisphere{ii}, '   ' titlEE{ii}])


        end
    case 5

  

        matDIR = dir('*.mat');
        matDIR2 = {matDIR.name};
        matDIR2 = transpose(matDIR2);

        timeLIn = matDIR2(contains(matDIR2,'TimeLine'));

        allDAyS = zeros(length(timeLIn),1);
        for ti = 1:length(timeLIn)

            load(timeLIn{ti},'outMAT')
            allDAyS(ti) = size(outMAT.month,2);

        end

        meanDAYS = mean(allDAyS); 
        stdDAYS = std(allDAyS);
        minDAYS = min(allDAyS);
        maxDAYS = max(allDAyS);

        array2table([meanDAYS,stdDAYS,minDAYS,maxDAYS],'VariableNames',{'Mean',...
            'STD','Min','Max'})

        table(timeLIn,allDAyS,'VariableNames',{'Subject_Hemi','Days'})

    case 6


        mainFig = figure;
        set(mainFig,'Position', [1485 311 1128 863]);
        tiledlayout(4,4,"Padding","tight");
        % Load data
        % GOOD SUBJECt
        % subjectID = '3';
        % hemisphere = 'L';


        % ARTIFACT 1
        % subjectID = '9';
        % hemisphere = 'L';

        % ARTIFACT 2
        % subjectID = '3';
        % hemisphere = 'L';

        % Non-ARTIFACT 1
        % subjectID = '8';
        % hemisphere = 'L';

        % Non-ARTIFACT 2
        % subjectID = '4';
        % hemisphere = 'R';

        subjectID = {'7','9','6','1'};
        hemisphere = {'L','L','R','R'};
        titlEE = {'8.8 Hz','8.8 Hz','12.7 Hz','10.7 Hz'};
        for ii = 1:4


            [tmData] = getPatDat(subjectID{ii} , hemisphere{ii} , 'TimeLine');

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            raw_LFP = unfurlLFP - (min(unfurlLFP));
            clean_LFP = raw_LFP;
            clean_LFP(clean_LFP > 2.2999e+09) = nan;
            % clean_LFP(clean_LFP < 50) = nan;
            mSunful = normalize(clean_LFP, 'range');
            smSunful = smoothdata(mSunful,'rloess',10,'omitnan');

            % Clean and prettified
            nexttile([1 4])
            cMAP = cividis;
            [reMAP] = reMapCmap(tmData,cMAP,smSunful,1,'timeBased');
            scatter(1:length(smSunful),smSunful,[],reMAP,'filled')

            maxVale = max(smSunful);
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
            title([subjectID{ii},'  ',hemisphere{ii}, '   ' titlEE{ii}])


        end


    case 7


        % mainFig = figure;
        % set(mainFig,'Position', [1485 311 1128 863]);
        % tiledlayout(4,4,"Padding","tight");

        % TABLE - MIN / MAX / Median / RAW
        % TABLE - MIN / MAX / Median / Normalized [range]

        subjectID = {'7','9','6','1'};
        hemisphere = {'L','L','R','R'};
        titlEE = {'8.8 Hz','8.8 Hz','12.7 Hz','10.7 Hz'};
        for ii = 1:4

            [tmData] = getPatDat(subjectID{ii} , hemisphere{ii} , 'TimeLine');

            nLFP = tmData.LFP;
            unfurlLFP = nLFP(:);
            raw_LFP = unfurlLFP - (min(unfurlLFP));
            clean_LFP = raw_LFP;
            clean_LFP(clean_LFP > 2.2999e+09) = nan;
            % clean_LFP(clean_LFP < 50) = nan;
            mSunful = normalize(clean_LFP, 'range');
            % smSunful = smoothdata(mSunful,'rloess',10,'omitnan');

            % % Clean and prettified
            % nexttile([1 4])
            % cMAP = cividis;
            % [reMAP] = reMapCmap(tmData,cMAP,smSunful,1,'timeBased');
            % scatter(1:length(smSunful),smSunful,[],reMAP,'filled')
            % 
            % maxVale = max(smSunful);
            % ylim([0 round(maxVale + 0.1,1)])
            % yticks([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
            % yticklabels([0 round((maxVale + 0.1)/2,2) round(maxVale + 0.1,1)])
            % ylabel('Scaled power')
            % dayStarts = round(linspace(1,length(smSunful)-144,length(smSunful)/144));
            % xticks(dayStarts);
            % xticklabels(1:length(smSunful)/144)
            % xlim([1, length(smSunful)])
            % xlabel('Days of recording')
            % set(gca,'TickLength',[0 .001])
            % title([subjectID{ii},'  ',hemisphere{ii}, '   ' titlEE{ii}])


        end

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
    if sum(matLog) == 0
        tmpDaOUT = nan;
    else
        load(matNames{matLog},'outEVENTS');
        tmpDaOUT = outEVENTS;
    end
elseif matches(tID, 'ProphetFOR')
    matDir = dir('*.csv');
    matNames = {matDir.name};
    matEls = split(matNames,'_');

    cIDs1 = matEls(:,:,1);
    cIDs2 = extractAfter(cIDs1,4);

    hIDs = matEls(:,:,2);

    tIDs1 = matEls(:,:,3);
    tIDs2 = extractBefore(tIDs1,'.');

    matLog = matches(cIDs2,cID) & matches(hIDs,hID) & matches(tIDs2,tID);
    outMAT = readtable(matNames{matLog});
    tmpDaOUT = outMAT;
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