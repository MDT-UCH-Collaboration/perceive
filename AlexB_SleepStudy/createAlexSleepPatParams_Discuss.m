function [] = createAlexSleepPatParams_Discuss(saveDIR)

% SPPD 11 - two dates
if nargin == 0
    saveDIR = 'E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';
end

patParams = struct;


%% 11 A
patParams.P11.A.L.OverSat = true;
patParams.P11.A.L.json = 'Report_Json_Session_Report_20220329T143105.json';
patParams.P11.A.R.OverSat = true;
patParams.P11.A.R.json = 'Report_Json_Session_Report_20220329T143105.json';
%% 11 B
patParams.P11.B.L.OverSat = true;
patParams.P11.B.L.json = 'Report_Json_Session_Report_20220317T110442.json';
patParams.P11.B.R.OverSat = true;
patParams.P11.B.R.json = 'Report_Json_Session_Report_20220317T114335.json';


%% Final step

dateCreate = datestr(now);

cd(saveDIR)
save('PatParmsDiscuss.mat','patParams','dateCreate');





