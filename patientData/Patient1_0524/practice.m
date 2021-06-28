%QUESTIONS: Dr. T, if your script shows mag v freq, why are we doing 
%spectrum to view in frequency domain? Aren't we there? 
%ask about offset time, aka why is it 3pm? %who knows! 
%      it is six hours 
%why is this only occuring on stim right, body left? will it always be 
%directional? why? 
%      assess stim on specifiic hemi, without spillage from other side 
%      need clean signal, hard with two stim sources
%      compare original w novel based on beta peak

%pwr in db is similar is mV like other one, db is similar enough
%electrode is showing oscillating voltage, freq also oscillates 
%power is output from transform fast fourier transform
%tells us magnitude at each frequency, we care about beta, bc science

%% loads in json file 
cd('C:\Users\sydne\Documents\MATLAB\ThompsonLab\Patient1_0524')
%initDir = dir('*.json');
jsonFiles = 'Report_Json_Session_Report_Patient1_right_init.json';

js = jsondecode(fileread(jsonFiles)); 

%%  
% figure;
% for i = 1:6 
%     testdata = js.LfpMontageTimeDomain(i).TimeDomainData;
%     [p,f] = pspectrum(testdata, 250, 'FrequencyLimits', [0 100]);
%     plot(f,db(p,"power"));
%     hold on; 
% end
% 
% xlabel("Frequency (Hz)")
% ylabel("Power Spectrum (dB)")
% legend("1", "2", "3", "4", "5", "6")
%% compare specific trials (is it junk?)
figure; 
    testdata = js.LfpMontageTimeDomain(1).TimeDomainData
    
    [p,f] = pspectrum(testdata);
    semilogx(f,p);
    semilogx(f,db(p,"power"));
    hold on;
    
    testdata = js.LfpMontageTimeDomain(7).TimeDomainData;
    [p,f] = pspectrum(testdata);
    semilogx(f,p);
    semilogx(f,db(p,"power"))
    hold on; 
    
    testdata = js.LfpMontageTimeDomain(13).TimeDomainData;
    [p,f] = pspectrum(testdata);
    semilogx(f,p);
    semilogx(f,db(p,"power"));
    hold on; 

    testdata = js.LfpMontageTimeDomain(20).TimeDomainData;
    [p,f] = pspectrum(testdata);
    semilogx(f,p);
    semilogx(f,db(p,"power"));
    hold on; 
    
    
xlabel("Frequency (Hz)")
ylabel("Power Spectrum (dB)")
legend("1", "7", "13", "20")
%% using function to compare specific trials 
close;
nonAvgTrialCompare(leng, js, 1, 7, 13, 19) 
%% 
%Generates differences between trials
testtime = {js.LfpMontageTimeDomain.FirstPacketDateTime};
cleantime = {}; %initialize new cell array for clean times
leng = numel(testtime); %have to use numel here bc length sucks? 
for i = 1:leng
    bettertime = replace(testtime{i}, {'T', 'Z'},  {' ', ''}); 
    cleantime{i} = datetime(bettertime) - hours(6); %converts string to datetime
    %fix by moving the time back six hours to make it correct 
end

%now find differentials between the times 
diff = {}; %q is this still a cell array?? 
for i = 1:(leng-1)
    diff{i} = cleantime{i+1} - cleantime{i};
end
%the first entry is the difference between cells 2 and 1, 
%second entry is diff between 3 and 2, etc 

%% run to get full graph
%{} gets the thing from the actual cell array
%() gets cell ARRAY! 

close
leng = numel({js.LfpMontageTimeDomain.Channel});
channels = unique({js.LfpMontageTimeDomain.Channel}); 
colors = 'krbgmc';
for c=1:length(channels)
    tempAvgPlot(7, leng, js, channels{c}, colors(c))
    hold on; 
end 

%add all the plot info 

xline(13)
xline(30)

%these are the order presented in channels

legendChan = {};
for i = 1:length(channels) 
    legendChan{i} = replace(channels{i}, '_', ' ');
end 

legend(legendChan{1}, legendChan{2}, legendChan{3}, legendChan{4}, legendChan{5}, legendChan{6}, "Beta lb", "Beta ub")
title("Freq vs. Power (db)")
xlabel("Frequency")
ylabel("Power") 

%% Direct compare two contact pairs 
close; 
channelA = channels{5}; 
channelB = channels{6};
directCompare(7, leng, js, channelA, channelB)

xline(13)
xline(30)
%these are the order presented in channels
legend(channelA, channelB, "Beta lb", "Beta ub")
title("Freq vs. Power (db)")
xlabel("Frequency")
ylabel("Power") 
