function [] = upV_dk_analysis(fileLOC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cd(fileLOC)
load('jattable_uPv.mat','mastertableJ')

% patient
patLIST = string(mastertableJ.patient);
patID = unique(patLIST);

allpatients = {};
allsides = {};
allcontacts = {};
allbetaPow = [];
allbetaFreq = [];
count = 1;
for pi = 1:length(patID)

    patTable = mastertableJ(matches(patLIST,patID(pi)),:);

    % side of brain
    sideID = unique(patTable.SideID);

    for si = 1:length(sideID)
        sideTable = patTable(matches(patTable.SideID,sideID(si)),:);

        conPairID = unique(sideTable.ChanID);

        for ci = 1:length(conPairID)

            conPtable = sideTable(matches(sideTable.ChanID,conPairID(ci)),:);

            maxBetaUpv = zeros(1,3);
            maxBetaHz = zeros(1,3);
            for fi = 1:height(conPtable)
                freqPow = conPtable.PF_Data{fi};
                betaLOC = freqPow.Frequency >= 13 & freqPow.Frequency <= 30;
                betaFrqs = freqPow.Frequency(betaLOC);
                [betaPow , betaPloc] = max(freqPow.Power(betaLOC));
                betaPHz = betaFrqs(betaPloc);
                maxBetaUpv(fi) = betaPow;
                maxBetaHz(fi) = betaPHz;
            end
            % Average max uPv

            % Average Hz location
            allpatients{count,1} = patID(pi);
            allsides{count,1} = sideID(si);
            allcontacts{count,1} = conPairID(ci);
            allbetaPow(count,1) = mean(maxBetaUpv);
            allbetaFreq(count,1) = mean(maxBetaHz);
            count = count + 1; 


        end
    end
end

finalTable = table(allpatients,allsides,allcontacts,allbetaPow,allbetaFreq,...
    'VariableNames',{'PatID','SideID','Hemi','BetaPw','BetaHz'})

end