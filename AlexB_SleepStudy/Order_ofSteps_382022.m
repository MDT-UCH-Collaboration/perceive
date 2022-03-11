%% Order of Scripts and Functions 
% -----Revise date - 3/8/2022
% Note: created

% GET ACT START AND STOP

% FIX 7 and 8 with ACT start and STOP


%% Step 1
% Process actigraphy
actigraphyProcess(subID)

%% Step 2
% Load in Pat struct
cd('D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020')
load("PatParms.mat","patParams");

%% Step 3 - Events
% One case
pat2use = 9;
side = 2; % for bilateral: 1 = 'L'

mainDIR = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
userDIRs = [mainDIR,'\Data\SPPD'];
userDIRe = '\JSON_LFP';
% saveDIRs = [mainDIR,'\Data\SPPD'];
saveDIRe = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
tabLOC = [mainDIR,'\summarydataTab.csv'];
actDloc = [mainDIR,'\Data\SPPD'];

patParsmsFs = fieldnames(patParams);
patParsmsFs = replace(patParsmsFs,'P','');
patFields = patParams.(['P',patParsmsFs{matches(patParsmsFs,num2str(pat2use))}]);
hemiFields = fieldnames(patFields);

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

%% Step 4 - Timeline
% One case
pat2use = 9;
side = 2; % for bilateral: 1 = 'L'

mainDIR = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
userDIRs = [mainDIR,'\Data\SPPD'];
userDIRe = '\JSON_LFP';
% saveDIRs = [mainDIR,'\Data\SPPD'];
saveDIRe = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
tabLOC = [mainDIR,'\summarydataTab.csv'];
actDloc = [mainDIR,'\Data\SPPD'];

patParsmsFs = fieldnames(patParams);
patParsmsFs = replace(patParsmsFs,'P','');
patFields = patParams.(['P',patParsmsFs{matches(patParsmsFs,num2str(pat2use))}]);
hemiFields = fieldnames(patFields);


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
    'saveDIR',saveDIR,'stagE',2,'userDIR',userDIR,...
    "actDIR",actDloc,'hemiS',hemiS,"tabLOC",tabLOC,"jsonDAT",jsoN)


%% Step 5 - Plot
close all
mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
sleepPlotHelp_fun(mainLOC , num2str(9))





