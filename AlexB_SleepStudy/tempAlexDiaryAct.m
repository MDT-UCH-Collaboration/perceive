% TO DO
% 4. ACTIGRAPHY - activity
% 5. ACTIGRAPHY - when 2 bed - when 2 rise


%% 

cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\SPPD1\DIARY_data')

tmpTab = readtable('SPPD_1_WT.xlsx');

intoBED = datestr(tmpTab.into_Bed(1));
lightsOUT = datestr(tmpTab.lights_out(1));
outOFbed = datestr(tmpTab.out_of_bed(1));


% What to do with it.

%%

cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\SPPD1\ACT_data')

tmpAct = readtable("1764467_SPPD1_4_8_2021_12_01_00_PM_New_Analysis.csv")



