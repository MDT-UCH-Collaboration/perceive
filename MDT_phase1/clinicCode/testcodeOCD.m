js = jsondecode(fileread('Report_Json_Session_Report_20220120T163554.json'));
data = js.BrainSenseTimeDomain.TimeDomainData; 
plot(data);
[p, f] = pspectrum(data); 
figure
plot(f, p);
figure
[~, ~, t] = pspectrum(data, 'spectrogram'); 
figure
waferfall(f, t, p')