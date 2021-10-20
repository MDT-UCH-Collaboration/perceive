%

load('SPPD3_TimeLine.mat')

nLFP = normalize(outMAT.LFP,'range');
% Plot raw
% plot(1:144,nLFP,'Color',[0.5 0.5 0.5])
mLFP = mean(nLFP,2);
sLFP = std(nLFP,[],2);
uLFP = mLFP + sLFP;
dLFP = mLFP - sLFP;

%%
hold on
% Plot mean
plot(1:144,mean(nLFP,2),'Color','r','LineWidth',2)
plot(1:144,uLFP,'Color',[1,0.5,0.5],'LineWidth',1,'LineStyle','--')
plot(1:144,dLFP,'Color',[1,0.5,0.5],'LineWidth',1,'LineStyle','--')


xl1 = xline(90,'-.','9 PM','DisplayName','Average Sales');
xl2 = xline(7,'-.','7 AM','DisplayName','Average Sales');
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'center';


xticks(1:144)
xticks(ceil(linspace(1,144,12)))
xticklabels(outMAT.TimeX(ceil(linspace(1,144,12))))
yticks([0 0.5 1])

ylabel('Normalized power - peak beta [32 Days]')
