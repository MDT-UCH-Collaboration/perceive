%% Order of Scripts and Functions 
% -----Revise date - 3/8/2022
% Note: created


%% Step 1
% Process actigraphy
actigraphyProcess(subID)

%% Step 2
% Load in Pat struct
cd('D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020')
load("PatParms.mat","patParams");

%% Step 3
% One case
pat2use = 10;

mainDIR = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
userDIRs = [mainDIR,'\Data\SPPD'];
userDIRe = '\JSON_LFP';
% saveDIRs = [mainDIR,'\Data\SPPD'];
saveDIRe = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
tabLOC = [mainDIR,'\summarydataTab.csv'];
actDloc = [mainDIR,'\Data\SPPD'];

patParsmsFs = fieldnames(patParams);
patFields = patParams.(patParsmsFs{contains(patParsmsFs,num2str(pat2use))});
hemiFields = fieldnames(patFields);

side = 2; % for bilateral: 1 = 'L'
% if length(hemiFields) == 1
patTAB = patFields.(hemiFields{side}); % for bilateral: 1 = 'L'
hemi = hemiFields{side}; % for bilateral: 1 = 'L'

hemiS = hemi;
patID = pat2use;
overSAT = patTAB.OverSat;
jsonNAMEs = patTAB.json;

saveDIR = saveDIRe;
userDIR = [userDIRs , num2str(patID) , userDIRe];
jsoN = jsonNAMEs;

perceive_sleepTandE_v6('overSAT',overSAT,'subID',patID,...
    'saveDIR',saveDIR,'stagE',1,'userDIR',userDIR,...
    "actDIR",actDloc,'hemiS',hemiS,"tabLOC",tabLOC,"jsonDAT",jsoN)