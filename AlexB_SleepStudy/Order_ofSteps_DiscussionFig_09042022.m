%% Order of Scripts and Functions 
% -----Revise date - 3/8/2022
% Note: created

% SPPD 11 - two dates


%% Step 1
% Process actigraphy
actigraphyProcess(subID,1)

%% Step 2
% Load in Pat struct
cd('E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020')
load("PatParmsDiscuss.mat","patParams");


%% Step 4 - Timeline
% One case
pat2use = 11;
side = 1; % for bilateral: 1 = 'L'
numT = 'A';

% mainDIR = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
mainDIR = 'E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
userDIRs = [mainDIR,'\DisucssionFigure\SPPD',filesep,numT];
userDIRe = '\JSON_LFP';
% saveDIRs = [mainDIR,'\Data\SPPD'];
% saveDIRe = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
saveDIRe = 'E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\DisucssionFigure\finaltest';
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

%% Step 5 - Run PyActrigraphy
% 1. Copy paste - raw CSV from Actiwatch to
% 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\ACTrawPy'
% 2. Execute python scrip - 'pyActOut.ipynb'


%% Step 6 - Process PyActigraphy
ml = 'C:\Users\Admin\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\';
sl = 'C:\Users\Admin\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
processACTpyData(ml,sl)

%% Step U - something weird happened to the first 2 days of Pat#11


%% Step 5 - Plot
close all
mainLOC = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\testSav';
sleepPlotHelp_fun(mainLOC , num2str(9))





