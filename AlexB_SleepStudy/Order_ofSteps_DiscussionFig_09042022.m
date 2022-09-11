


%% Order of Scripts and Functions 
% -----Revise date - 3/8/2022
% Note: created

% SPPD 11 - two dates


%% Step 1
% Process actigraphy
actigraphyProcess(subID,1)

%% Step 2
% Load in Pat struct
cd('C:\Users\Admin\Desktop\sppd11tesat')
firstj = 'Report_Json_Session_Report_20220317T114335.json';
jsonfir = jsondecode(fileread(firstj));

secondj = 'Report_Json_Session_Report_20220329T143105.json';
jsonsec = jsondecode(fileread(secondj));


%% Step 4 - Timeline
% SPPD11 L First session
pat2use = 11;
side = 'R';
session = '2';

mainDIR = 'E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
saveDIR = 'C:\Users\Admin\Desktop\sppd11tesat\saveOUT';
userDIR = 'C:\Users\Admin\Desktop\sppd11tesat';
tabLOC = [mainDIR,'\summarydataTab.csv'];

hemiS = side;
patID = pat2use;
jsoN = secondj;

perceive_sleepTandE_Discussion('overSAT',1,'subID',patID,...
    'saveDIR',saveDIR,'stagE',2,'userDIR',userDIR,...
    'hemiS',hemiS,"tabLOC",tabLOC,"jsonDAT",jsoN,"sessNUM",session)

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





