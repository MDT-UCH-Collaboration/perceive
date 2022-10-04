function [] = mdt_p1_p5figure(dbdir,patid)

% C:\Users\John\Documents\GitHub\perceive\MDT_Phase1_patient5
cd(dbdir)


if patid == 5
    load("patient5_LeftBSS.mat","outTABLE")
elseif patid == 4
    load("patient4_RightBSS.mat","outTABLE")
end

uniChan = unique(outTABLE.ChanID);

allmeans = zeros(height(outTABLE.PF_Data{1}),size(uniChan,1));
chanID = {'1-3','1-2','2-3','0-1','0-3','0-2'};
for ci = 1:size(uniChan,1)

    tmpChan = uniChan{ci};
    chanInd = matches(outTABLE.ChanID,tmpChan);
    chanTab = outTABLE(chanInd,:);

    allsess = zeros(height(chanTab.PF_Data{1}),3);
    for ri = 1:3
        allsess(:,ri) = chanTab.PF_Data{ri}.Power;
    end
    meanSess = mean(allsess,2);
    allmeans(:,ci) = meanSess;

end

freqXaxis = outTABLE.PF_Data{1}.Frequency;

% smooth data
allmeanSM = smoothdata(allmeans,'gaussian',150);
% Trim to 60
allmeanSM60 = allmeanSM(freqXaxis <= 60,:);
freqAxis60 = freqXaxis(freqXaxis <= 60);

% Reorder legend by max beta
albetaP = zeros(6,1);
for bi = 1:6
    tmpAsm = allmeanSM60(:,bi);
    tmpBeta = tmpAsm(freqAxis60 >= 13 & freqAxis60 <= 30);
    albetaP(bi) = max(tmpBeta);
end
[~,maxBsort] = sort(albetaP);

colors_p = [121, 125, 98;
           255, 243, 176;
           241, 220, 167;
           255, 203, 105;
          208, 140, 96;
          153, 123, 102];
colors_pI = colors_p/255;

for pi = 1:6
    plot(freqAxis60,allmeanSM60(:,maxBsort(pi)),"Color",colors_pI(pi,:),'LineWidth',2)
    hold on
end


% Plot xline for beta and theta
xline([13 30],'-',{'beta',' '})
xline([4 8],'--',{'theta',' '})

legend(chanID(maxBsort))

xticks([0 30 60])
xlabel('Frequency (Hz)')

ylim([0 1.15])
yticks([0 0.5 1 1.15])
ylabel('Normalized uPv')

axis square





end