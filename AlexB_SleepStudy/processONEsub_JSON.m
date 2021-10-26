
userDIRs = 'C:\Users\Admin\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\SPPD';
userDIRe = '\JSON_LFP';
saveDIRs = 'C:\Users\Admin\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\SPPD';
saveDIRe = '\MAT_data';
tabLOC = "C:\Users\Admin\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\summarydataTab.csv";

jsonNAMEs = {'Report_Json_Session_Report_20210521T115931[1].json',...
             'Report_Json_Session_Report_20210521T105436[1].json',...
             'Report_Json_Session_Report_20210514T110817.json',...
             'Report_Json_Session_Report_20210615T131452.json',...
             'Report_Json_Session_Report_20210707T111453.json',...
             'Report_Json_Session_Report_20210816T120043.json',...
             'Report_Json_Session_Report_20210726T094956.json',...
             'Report_Json_Session_Report_20210726T100116.json'};


%%
hemiS = {'R','L','L','R','R','R','L','R'};
patID = [1, 2, 3, 4, 5, 6, 7, 7];
overSAT = [true , false , true , true , true, true, true, true];
 
for pi = 1:8

    overSATi = overSAT(pi);
    saveDIR = [saveDIRs , num2str(patID(pi)) , saveDIRe];
    userDIR = [userDIRs , num2str(patID(pi)) , userDIRe];
    jsoN = jsonNAMEs{pi};

    perceive_sleepTandE_v2('overSAT',overSATi,'subID',patID(pi),...
        'saveDIR',saveDIR,'stagE',2,'userDIR',actDIR,...
        userDIR,'hemiS',hemiS(pi),"tabLOC",tabLOC,"jsonDAT",jsoN)

end