%
function [] = sleepPlotHelp_fun(mainLOC)
% C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\SummaryMAT

cd(mainLOC)

matDir = dir('*.mat');
matDir2 = {matDir.name};

for i = 1:length(matDir)

    load(matDir2{i},'outMAT')

    nLFP = normalize(outMAT.LFP,'range');
    % Plot raw
    % plot(1:144,nLFP,'Color',[0.5 0.5 0.5])
    mLFP = mean(nLFP,2);
    sLFP = std(nLFP,[],2);
    uLFP = mLFP + sLFP;
    dLFP = mLFP - sLFP;

    hold on
    % Plot mean
    plot(1:144,mean(nLFP,2),'Color','r','LineWidth',2)
    plot(1:144,uLFP,'Color',[1,0.5,0.5],'LineWidth',1,'LineStyle','--')
    plot(1:144,dLFP,'Color',[1,0.5,0.5],'LineWidth',1,'LineStyle','--')

    xl1 = xline(90,'-.','9 PM','DisplayName','Average Sales');
    xl2 = xline(7,'-.','7 AM','DisplayName','Average Sales');
    xl1.LabelVerticalAlignment = 'top';
    xl1.LabelHorizontalAlignment = 'center';
    xl2.LabelVerticalAlignment = 'top';
    xl2.LabelHorizontalAlignment = 'center';

    xticks(1:144)
    xticks(ceil(linspace(1,144,12)))
    xticklabels(outMAT.TimeX(ceil(linspace(1,144,12))))
    yticks([0 0.5 1])

    daYS = size(nLFP,2);

    ylabel(['Normalized power - peak beta [', num2str(daYS) ,' Days]'])

    pause
    close all

end