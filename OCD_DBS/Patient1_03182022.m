%% Read JSON File
cd('C:\Users\johna\Desktop\ritalin_LFP\PRE')
subjectFILE = 'Report_Json_Session_Report_20220317T110442.json';
js = jsondecode(fileread(subjectFILE));
%%
streamF = js.BrainSenseTimeDomain(1).TimeDomainData;

%% highpass
streamH = highpass(streamF,0.5);
%%
streamN = spectrumInterpolation(streamH, 250, 60, 4, 2);

%% 
plot(streamN)

