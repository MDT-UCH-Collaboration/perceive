cd('C:\Users\sydne\Documents\github\perceive\patientData\Patient3_0630')
%jsonFiles = 'Report_Json_Session_Report_20210630T155026.json';
%file above is left side brainsense survey

%left brain streaming json ending in 3636 
%this file spcifically is 3398 kb
jsonFiles = 'Report_Json_Session_Report_20210630T163636.json';
js1 = jsondecode(fileread(jsonFiles));
%diagnostic even data for 6/30 starts at index 18 

%left brain streaming json ending in 3648 
%this file is 12 kb :( (sad and gross) 
%this file lacks any brainsense data 
jsonFiles = 'Report_Json_Session_Report_20210630T163648.json';
js2 = jsondecode(fileread(jsonFiles));
%diagnostic event data for 6/30 starts at index 18 
%diff between session start and end is 15 seconds 

%sets it back to the path where all the other functions are
cd('C:\Users\sydne\Documents\github\perceive\postAnalysis')

%note: highest beta is 
%js1_ma = js1.Groups.Initial.ProgramSettings.LeftHemisphere.Programs.AmplitudeInMilliAmps; 
%js1_pw = js1.Groups.Initial.ProgramSettings.LeftHemisphere.Programs.PulseWidthInMicroSecond;
%I think this is the streaming info for the entry program
%?? bc the active group is true 
%stim field tells us that stim is off?? 

%from the white paper from Michelle 
%brainsensetimedomain is time domain data from when streaming button is
%pressed 
%brainsenselfp is power data from when streaming is pressed 

%%
%this code makes a struct that I can selectively store the relevant info of
%from the json streaming file 

data = struct; 
data.t1.channel = js1.BrainSenseLfp(1).Channel; 
data.t1.time = datetime(replace({js1.BrainSenseLfp(1).FirstPacketDateTime}, {'T', 'Z'},  {' ', ''}));
data.t1.pw = js1.BrainSenseLfp(1).TherapySnapshot.Left.PulseWidthInMicroSecond; 
data.t1.ma = js1.BrainSenseLfp(1).TherapySnapshot.Left.LowerLimitInMilliAmps; 
data.t1.timeDomain = js1.BrainSenseTimeDomain(1).TimeDomainData; 

%gather the werido lfp data 
%data.t1.lfp = js1.BrainSenseLfp(1).LfpData(1).Left.LFP; 
lfp_leng = numel(js1.BrainSenseLfp(1).LfpData); 
lfp = zeros(lfp_leng, 1);
for i = 1:lfp_leng 
    lfp(i) = js1.BrainSenseLfp(1).LfpData(i).Left.LFP; 
end 
data.t1.lfp = lfp;

[p,f] = pspectrum(data.t1.timeDomain, 250, 'FrequencyLimits', [0 100]); 
plot(f, p); 
title("freq vs pwr"); 
xlabel("freq (Hz)");
ylabel("pwr"); 

spec = spectrogram(data.t1.timeDomain); 
plot(spec)
%% 
%making time lmoa
figure; 
ticks = js1.BrainSenseTimeDomain(1).TicksInMses; 
ticks = str2num(ticks);  
ticks = ticks';  

packs = str2num(js1.BrainSenseTimeDomain(1).GlobalPacketSizes); 
pack_size = packs(1); 

smpl_rate = js1.BrainSenseTimeDomain.SampleRateInHz; 

timeleng = length(ticks);  

TDtime = ticks(timeleng) - ((pack_size-1)/smpl_rate):1/smpl_rate:ticks(timeleng); 
for i = 1:timeleng-1
    prev_pack = ticks(timeleng-i) - ((pack_size-1)/smpl_rate):1/smpl_rate:ticks(timeleng-i); 
    TDtime = [prev_pack, TDtime];  
end 

TDtime = TDtime - (ticks(1) - ((pack_size-1)/smpl_rate)); 
TDtime = TDtime/1000; 

%TDtime is the appropriate time in seconds!!!!!!!!!!!!!!!!! :) I think! 

data.t1.time = TDtime'; 

%plot(data.t1.time, data.t1.timeDomain, 'Marker', 'x', 'LineStyle', 'none')

%% 
% %spectrogram? 
% y = diff(data.t1.time);
% yy = diff(y);
% pspectrum(data.t1.timeDomain, data.t1.time, 'spectrogram')