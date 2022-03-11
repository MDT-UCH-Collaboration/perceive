function [] = createAlexSleepPatParams(saveDIR)

if nargin == 0
    saveDIR = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
end

patParams = struct;

%% 1
patParams.P1.R.OverSat = true;
patParams.P1.R.json = 'Report_Json_Session_Report_20210521T115931[1].json';
%% 2
patParams.P2.L.OverSat = false;
patParams.P2.L.json = 'Report_Json_Session_Report_20210521T105436[1].json';
%% 3
patParams.P3.L.OverSat = true;
patParams.P3.L.json = 'Report_Json_Session_Report_20210514T110817.json';
%% 4
patParams.P4.R.OverSat = true;
patParams.P4.R.json = 'Report_Json_Session_Report_20210615T131452.json';
%% 5
patParams.P5.R.OverSat = true;
patParams.P5.R.json = 'Report_Json_Session_Report_20210707T111453.json';
%% 6
patParams.P6.R.OverSat = true;
patParams.P6.R.json = 'Report_Json_Session_Report_20210816T120043.json';
%% 7
patParams.P7.L.OverSat = true;
patParams.P7.L.json = 'Report_Json_Session_Report_20210726T094956.json';
patParams.P7.R.OverSat = true;
patParams.P7.R.json = 'Report_Json_Session_Report_20210726T100116.json';
%% 8
patParams.P8.L.OverSat = true;
patParams.P8.L.json = 'Report_Json_Session_Report_20211123T131149.json';
%% 9
patParams.P9.L.OverSat = true;
patParams.P9.L.json = 'Report_Json_Session_Report_20220104T121534.json';
patParams.P9.R.OverSat = true;
patParams.P9.R.json = 'Report_Json_Session_Report_20220104T121534.json';
%% 10
patParams.P10.L.OverSat = true;
patParams.P10.L.json = 'Report_Json_Session_Report_20220124T081628.json';
patParams.P10.R.OverSat = true;
patParams.P10.R.json = 'Report_Json_Session_Report_20220124T081944.json';

dateCreate = datestr(now);

cd(saveDIR)
save('PatParms.mat','patParams','dateCreate');





