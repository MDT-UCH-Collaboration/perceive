cd('C:\Users\Admin\Documents\Github\perceive\postAnalysis')

dataTAB = readtable("anovaDataAlljat3.xlsx");

patUNI = unique(dataTAB.Patient);

coloRS = [39, 84, 66;...
         40, 114, 113;...
         42, 157, 143;...
         138, 177, 125;...
         233, 196, 106;...
         239, 179, 102;...
         244, 162, 97;...
         238, 137, 89;...
         231, 111, 81];

coloRSrgb = coloRS/255;
colCOUNT = 1;
for pi = 1:length(patUNI)

    tmpPI = dataTAB(matches(dataTAB.Patient,patUNI(pi)),:);

    hemiUNI = unique(tmpPI.Hemi);

    for hi = 1:length(hemiUNI)

        tmpHE = tmpPI(matches(tmpPI.Hemi,hemiUNI(hi)),:);

        blankCanvas = zeros(height(tmpHE),3);

        % black non-sig
        nonSIGind = tmpHE.adj_p > 0.05;
        % gray sig no match
        sigNoMatch = tmpHE.adj_p < 0.05 & ~matches(tmpHE.MaxID,tmpHE.maxBcPair(1));
        blankCanvas(sigNoMatch,:) = repmat([0.5 0.5 0.5],sum(sigNoMatch),1);
        trueMAtch = tmpHE.adj_p < 0.05 & matches(tmpHE.MaxID,tmpHE.maxBcPair(1));
        blankCanvas(trueMAtch,:) = repmat(coloRSrgb(colCOUNT,:),sum(trueMAtch),1);

        meanNS = mean(tmpHE.adj_p(nonSIGind));
        stdNS = std(tmpHE.adj_p(nonSIGind));
        upstd = meanNS + stdNS;
        dnstd = meanNS - stdNS;

        mX = zeros(1,1) + colCOUNT;
        scatter(mX,meanNS,70,'black','filled')
        line([mX mX],[dnstd upstd],'Color','k','LineWidth',1)
        hold on

        canVas2use = blankCanvas(:,1) ~= 0;
        canVas2use2 = blankCanvas(canVas2use,:);
        pAdj = tmpHE.adj_p(canVas2use);

        xVals = zeros(height(canVas2use2),1) + colCOUNT;
        s = scatter(xVals,pAdj,40,canVas2use2,'filled');
        s.XJitter = "density";
        s.XJitterWidth = 0.3;

        colCOUNT = colCOUNT + 1;

    end
end
xlim([0 10])
ylim([0 .2])





% Cross corr with max beta power

coloRSrgb = coloRS/255;
colCOUNT = 6;
for pi = 1:length(patUNI)

    tmpPI = dataTAB(matches(dataTAB.Patient,patUNI(pi)),:);

    hemiUNI = unique(tmpPI.Hemi);

    for hi = 1:length(hemiUNI)

        tmpHE = tmpPI(matches(tmpPI.Hemi,hemiUNI(hi)),:);

        blankCanvas = zeros(height(tmpHE),3);

        % black non-sig
        nonSIGind = tmpHE.adj_p > 0.05;
        % gray sig no match
        sigNoMatch = tmpHE.adj_p < 0.05 & ~matches(tmpHE.MaxID,tmpHE.maxBcPair(1));
        blankCanvas(sigNoMatch,:) = repmat([0.5 0.5 0.5],sum(sigNoMatch),1);
        trueMAtch = tmpHE.adj_p < 0.05 & matches(tmpHE.MaxID,tmpHE.maxBcPair(1));
        blankCanvas(trueMAtch,:) = repmat(coloRSrgb(colCOUNT,:),sum(trueMAtch),1);

        NSall = tmpHE.adj_p(nonSIGind);
        MBall = zeros(height(tmpHE),1);
        for mbi = 1:height(tmpHE)
            tmpID = tmpHE.MaxID{mbi};
            matLOC = find(matches(table2array(tmpHE(mbi,{'GroupA','GroupB'})),tmpID));
            switch matLOC
                case 1
                      MBall(mbi) = tmpHE.MaxB_GA(mbi);
                case 2
                      MBall(mbi) = tmpHE.MaxB_GB(mbi);
            end
        end

        MBnsAll = MBall(nonSIGind);
        scatter(NSall,MBnsAll,70,'black','filled')
   
        hold on

        canVas2use = blankCanvas(:,1) ~= 0;
        canVas2use2 = blankCanvas(canVas2use,:);
        pAdj = tmpHE.adj_p(canVas2use);
        MBsig = MBall(canVas2use);

        s = scatter(pAdj,MBsig,70,canVas2use2,'filled');
        s.XJitter = "density";
        s.XJitterWidth = 0.3;

    end
end
xlabel('FDR corrected p-values')
ylabel('Max Beta power')
xlim([0 10])
ylim([0 .2])





% Cross corr with max beta power

coloRSrgb = coloRS/255;
colCOUNT = 1;
allfracs = zeros(9,2);
legendID = cell(9,1);
for pi = 1:length(patUNI)

    tmpPI = dataTAB(matches(dataTAB.Patient,patUNI(pi)),:);

    hemiUNI = unique(tmpPI.Hemi);

    for hi = 1:length(hemiUNI)

        legPI = ['S',patUNI{pi}, ' ' , hemiUNI{hi}(1)];
        legendID{colCOUNT} = legPI;

        tmpHE = tmpPI(matches(tmpPI.Hemi,hemiUNI(hi)),:);

        blankCanvas = zeros(height(tmpHE),3);

        % black non-sig
        nonSIGind = tmpHE.adj_p > 0.05;
        % gray sig no match
        sigNoMatch = tmpHE.adj_p < 0.05 & ~matches(tmpHE.MaxID,tmpHE.maxBcPair(1));
        blankCanvas(sigNoMatch,:) = repmat([0.5 0.5 0.5],sum(sigNoMatch),1);
        trueMAtch = tmpHE.adj_p < 0.05 & matches(tmpHE.MaxID,tmpHE.maxBcPair(1));
        blankCanvas(trueMAtch,:) = repmat(coloRSrgb(colCOUNT,:),sum(trueMAtch),1);

%         fracNS = sum(nonSIGind)/height(tmpHE);
        fracSigNotMB = sum(sigNoMatch)/height(tmpHE);
        fracSigMB = sum(trueMAtch)/height(tmpHE);

        allfracs(colCOUNT,:) = [fracSigNotMB , fracSigMB];

%         mX = zeros(1,1) + colCOUNT;
%         scatter(mX,fracNS,70,'black','filled')
%         hold on
%         scatter(mX + 0.1,fracSigNotMB,70,[0.5 0.5 0.5],'filled')
%         s = scatter(mX + 0.2,fracSigMB,70,coloRSrgb(7,:),'filled');
        colCOUNT = colCOUNT + 1;

    end
end

bar(allfracs)
xlabel('Subject / hemipshere')
ylabel('Fraction of pair-wise comparisons')
legend({'Sig. All other contacts','Sig. MB contact'})
xticklabels(legendID)
axis square
