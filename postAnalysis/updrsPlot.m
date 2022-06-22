cd('D:\Dropbox\MDT_PhaseI_Manuscripts')

upTab = readtable('updrsSum.xlsx');

offD = matches(upTab.cond,'OFF');
entD = matches(upTab.cond,'Entry');
expD = matches(upTab.cond,'Exp');
expDa = matches(upTab.cond,'Exp') & matches(upTab.ExpType, 'maxb');
expDb = matches(upTab.cond,'Exp') & ~matches(upTab.ExpType, 'maxb');

line([zeros(1,sum(offD)) ; ones(1,sum(offD))] ,...
    [transpose(upTab.TOTAL(offD)) ; transpose(upTab.TOTAL(entD))],...
    'Color',[100/255 0/255 162/255],'LineWidth',1.5);
hold on
line([ones(1,sum(offD)) ; ones(1,sum(offD))+1] ,...
    [transpose(upTab.TOTAL(entD)) ; transpose(upTab.TOTAL(expD))],...
    'Color',[0/255 13/255 115/255],'LineWidth',1.5);

s1 = scatter(zeros(sum(offD),1),upTab.TOTAL(offD),90,'r','filled');
s1.MarkerFaceAlpha = 0.9;
s2 = scatter(ones(sum(offD),1),upTab.TOTAL(entD),90,'b','filled');
s2.MarkerFaceAlpha = 0.9;


% s3 = scatter(ones(sum(offD),1)+1,upTab.TOTAL(expD),90,'k','filled');
s3a = scatter(ones(sum(expDa),1)+1,upTab.TOTAL(expDa),90,'k','filled');
s3b = scatter(ones(sum(expDb),1)+1,upTab.TOTAL(expDb),90,'k');

% s3 = scatter(ones(sum(offD),1)+1,upTab.TOTAL(expD),90,'k','filled');
s3a.MarkerFaceAlpha = 0.9;

xticks([0 1 2])
xticklabels({'Off','Entry','Exp.'})
xlim([-0.2 2.2])
yticks([0 13 26])
ylim([0 26])
ylabel('MDS-UPDRS Total Score')


axis square

% omni
text(1,25, 'F = 7.9, (2,24), p = 0.002')
% off v entry
text(-0.1,23, 'oVen 3.0, 0.01')
% off v exp
text(0.7,23, 'oVex 3.7, 0.003')
% entry v exp
text(1.6,23, 'enVex 0.66, 0.78')

