clear variables;clc; close all

% Import and read file
[val,path] = Import_Read;

% Generate/Read variables
fs=250;
nfft=250;
window=250;
overlap=150;

sensechannelcolumnsindex = {'TipOffsetDef.ZERO',1;
    'TipOffsetDef.FOUR',2;
    'TipOffsetDef.EIGHT',3;
    'TipOffsetDef.TWELVE',4};

channelelectrodemap = {
    'ZERO_THREE';
    'ONE_THREE';
    'ZERO_TWO';
    'ZERO_AND_ONE';
    'ONE_AND_TWO';
    'TWO_AND_THREE'};

sensechannelmapped = {
    "0-3","4-7","8-11","12-15";
    "1-3","5-7","9-11","13-15";
    "0-2","4-6","8-10","12-14";
    "0-1","4-5","8-9","12-13";
    "1-2","5-6","9-10","13-14";
    "2-3","6-7","10-11","14-15"};

hemierasestring = 'HemisphereLocationDef.';
try
    electrodemap = {erase(string(val.LeadConfiguration.Final(1).Hemisphere),hemierasestring),string(val.LeadConfiguration.Final(1).TipOffset);
        erase(string(val.LeadConfiguration.Final(2).Hemisphere),hemierasestring),string(val.LeadConfiguration.Final(2).TipOffset)};
catch
    electrodemap = {erase(string(val.LeadConfiguration.Final(1).Hemisphere),hemierasestring),string(val.LeadConfiguration.Final(1).TipOffset)};
end

leftchannels = [];
rightchannels = [];


figure('name','BrainSenseSurvey (Pwelch)')

for i = 1:min([length(val.LfpMontageTimeDomain),12]) % for each stream i
    for j = 1:size(electrodemap,1) % for each hemisphere j
        if contains(string(val.LfpMontageTimeDomain(i).Channel),electrodemap{j,1},'IgnoreCase',true) % identify hemisphere
            for k = 1:4 % for each tip offset
                if electrodemap{j,2} == sensechannelcolumnsindex{k,1}
                    sensechannelcolumn = sensechannelcolumnsindex{k,2};
                end
            end
            for k = 1:6 % for each electrode pair
                if contains(string(val.LfpMontageTimeDomain(i).Channel),channelelectrodemap{k})
                    electrodepair = sensechannelmapped{k,sensechannelcolumn};
                    if contains(string(val.LfpMontageTimeDomain(i).Channel),'LEFT')
                        leftchannels = [leftchannels,electrodepair];
                        subplot(2,1,1)
                        
                        data=val.LfpMontageTimeDomain(i).TimeDomainData;
                        [pxx,f] = pwelch(data,window,overlap,nfft,fs);
                        semilogy(f,sqrt(pxx))
                        hold on
                    else
                        rightchannels = [rightchannels,electrodepair];
                        subplot(2,1,2)
                        
                        data=val.LfpMontageTimeDomain(i).TimeDomainData;
                        [pxx,f] = pwelch(data,window,overlap,nfft,fs);
                        semilogy(f,sqrt(pxx))
                        hold on
                    end
                    break
                end
            end
        end
    end
    subplot(2,1,1)
    ylabel('PSD (uV/rthz)')
    title('BrainSense Survey (Left)')
    ylim([.01 10])
    legend(leftchannels)
    
    subplot(2,1,2)
    xlabel('Frequency (Hz)')
    ylabel('PSD (uV/rthz)')
    title('BrainSense Survey (Right)')
    ylim([.01 10])
    legend(rightchannels)
    
    %     contains(string(val.LfpMontageTimeDomain(i).Channel),channelelectrodemap);
end
