function perceive_vMDTclinictest(inPS)
% https://github.com/neuromodulation/perceive
% v0.1 Contributors Wolf-Julian Neumann, Tomas Sieger, Gerd Tinkhauser
% This is an open research tool that is not intended for clinical purposes.

arguments
    inPS.userPC (1,1) string = "JATwork"
    inPS.subID (1,:) char = '001' % IN USE for CASE NUMBER
    inPS.postN (1,1) double = 1
    inPS.userDIR (1,1) string = "NA"
    inPS.seluDIR (1,1) logical = true
    inPS.saveDIR (1,1) string = "NA"
    inPS.selsDIR (1,1) logical = true
    inPS.anlyz (1,1) double = 1
    inPS.studY (1,:) char = '20-2508'
    inPS.pltCl (1,1) logical = 0
end


%% OUTPUT
% The script generates BIDS inspired subject and session folders with the
% ieeg format specifier. All time series data are being exported as
% FieldTrip .mat files, as these require no additional dependencies for creation.
% You can reformat with FieldTrip and SPM to MNE
% python and other formats (e.g. using fieldtrip2fiff([fullname '.fif'],data))

%% Recording type output naming
% Each of the FieldTrip data files correspond to a specific aspect of the
% Recording session:
% LMTD = LFP Montage Time Domain - BrainSenseSurvey
% IS = Indefinite Streaming - BrainSenseStreaming
% CT = Calibration Testing - Calibration Tests
% BSL = BrainSense LFP (2 Hz power average + stimulation settings)
% BSTD = BrainSense Time Domain (250 Hz raw data corresponding to the BSL
% file)

%% TODO:
% ADD DEIDENTIFICATION OF COPIED JSON
% BUG FIX UTC?
% ADD BATTERY DRAIN
% ADD BSL data to BSTD ephys file
% ADD PATIENT SNAPSHOT EVENT READINGS
% IMPROVE CHRONIC DIAGNOSTIC READINGS
% ADD Lead DBS Integration for electrode location

if inPS.seluDIR && strcmp(inPS.userDIR,"NA")
    [fileDIR] = uigetdir();
else
    fileDIR = inPS.userDIR;
end

if inPS.selsDIR && strcmp(inPS.saveDIR,"NA")
    [saveLOC] = uigetdir();
else
    saveLOC = inPS.saveDIR;
end

cd(fileDIR)
initDir = dir('*.json');
jsonFiles = {initDir.name};

for jsoni = 1:length(jsonFiles)
    
    filename = jsonFiles{jsoni};
    disp(['RUNNING ' filename])
    
    js = jsondecode(fileread(filename));
    % CLEAR PATIENT INFO - DEIDENTIFY
    disp('Clear patient info')
    try
        js.PatientInformation.Initial.PatientFirstName ='';
        js.PatientInformation.Initial.PatientLastName ='';
        js.PatientInformation.Initial.PatientDateOfBirth ='';
        js.PatientInformation.Final.PatientFirstName ='';
        js.PatientInformation.Final.PatientLastName ='';
        js.PatientInformation.Final.PatientDateOfBirth ='';
    catch
        js = rmfield(js,'PatientInformation');
        js.PatientInformation.Initial.PatientFirstName ='';
        js.PatientInformation.Initial.PatientLastName ='';
        js.PatientInformation.Initial.PatientDateOfBirth ='';
        js.PatientInformation.Initial.Diagnosis ='';
        js.PatientInformation.Final.PatientFirstName ='';
        js.PatientInformation.Final.PatientLastName ='';
        js.PatientInformation.Final.PatientDateOfBirth ='';
        js.PatientInformation.Final.Diagnosis = '';
    end
    
    % Initial data to extract
    
    infofields = {'SessionDate','SessionEndDate','PatientInformation',...
        'DeviceInformation','BatteryInformation','LeadConfiguration',...
        'Stimulation','Groups','Stimulation','Impedance','PatientEvents',...
        'EventSummary','DiagnosticData','LFPMontage'};
    
    for availINfld = 1:length(infofields)
        if isfield(js,infofields{availINfld})
            hdr.(infofields{availINfld})=js.(infofields{availINfld});
        end
    end
    
    disp(js.DeviceInformation.Final.NeurostimulatorLocation)
    
    [hdr.SessionDate] = fix2time(js.SessionDate , 1);
    
    if ~isempty(js.PatientInformation.Final.Diagnosis)
        hdr.Diagnosis = strsplit(js.PatientInformation.Final.Diagnosis,'.');
        hdr.Diagnosis = hdr.Diagnosis{2};
        if strcmp(hdr.Diagnosis,'ParkinsonsDisease')
            hdr.Diagnosis = 'PD';
        end
    else
        hdr.Diagnosis = '';
    end
    
    hdr.OriginalFile = filename;
    
    [hdr.ImplantDate] = fix2time(js.DeviceInformation.Final.ImplantDate , 1);
    hdr.BatteryPercentage = js.BatteryInformation.BatteryPercentage;
    leadLOC = strsplit(hdr.LeadConfiguration.Final(1).LeadLocation,'.');
    hdr.LeadLocation = leadLOC{2};
    hemiLOC = strsplit(hdr.LeadConfiguration.Final(1).Hemisphere,'.');
    hdr.Hemisphere = hemiLOC{2};
    
    hdr.subject = ['sub-',inPS.subID,'_', inPS.studY,'_', hdr.LeadLocation, '_', hdr.Hemisphere ];
    
    hdr.session = ['ses-' char(datetime(hdr.SessionDate,'format','yyyyMMddhhmmss')) num2str(hdr.BatteryPercentage)];
    
    if ~exist(fullfile(saveLOC,hdr.subject,hdr.session,'ieeg'),'dir')
        mkdir(fullfile(saveLOC,hdr.subject,hdr.session,'ieeg'));
    end
    hdr.fpath = fullfile(saveLOC,hdr.subject,hdr.session,'ieeg');
    hdr.fname = [hdr.subject '_' hdr.session];
    hdr.chan = ['LFP_' hdr.LeadLocation];
    hdr.d0 = datetime(js.SessionDate(1:10));
    hdr.js = js;
    if ~exist('datafields','var')
        datafields = sort({'EventSummary','Impedance','MostRecentInSessionSignalCheck',...
            'BrainSenseLfp','BrainSenseTimeDomain','LfpMontageTimeDomain',...
            'IndefiniteStreaming','BrainSenseSurvey','CalibrationTests','PatientEvents',...
            'DiagnosticData','LFPMontage'});
    end
    
    useDF = 'LFPMontage';
    switch useDF
        case 'LFPMontage'
            data = js.(useDF);
            
            cdata = {data(:).SensingElectrodes}';
            cdata2 = cellfun(@(x) strsplit(x,'.'), cdata,'UniformOutput',false);
            cdata3 = cellfun(@(x) x{2}, cdata2,'UniformOutput',false);
            extData = replace(cdata3,{'ZERO','ONE','TWO','THREE','_AND','_'},{'0','1','2','3','',''});
            
            peakMV = [data.PeakMagnitudeInMicroVolt];
            [~, maxIND] = max(peakMV); 
            figure;
            for i = 1:6
                if i == maxIND
                    plot(data(i).LFPFrequency ,data(i).LFPMagnitude,'LineWidth',3)
                else
                plot(data(i).LFPFrequency ,data(i).LFPMagnitude,'LineWidth',1)
                end
                hold on
            end
            
            
            xline(13,'--')
            xline(30,'--')
            legend(extData)

            
    end
end



hdr.DeviceInformation.Final.NeurostimulatorLocation
end







function [timeFix] = fix2time(inTIME , procStep)

switch procStep
    case 1
        timeFix = datetime(replace(inTIME(1:end-1),{'T'},{' '}));
        
    case 2
        timeFix = replace(inTIME(1:end-1),{'T'},{' '});
        
end

end



function [lfpsettings , stimchannels] = hemiCheck(inDAT,procStep)

switch procStep
    case 1
        
        rls = {'Left','Right'};
        stimName = {'STIM_L_','STIM_R_'};
        lfpsettings = cell(2,1);
        stimchannels = cell(1,2);
        for rl = 1:2
            
            if isfield(inDAT,rls{rl})
                tmpChanIDs = inDAT.(rls{rl});
                lfpsettings{rl,1} = ['PEAK' num2str(round(tmpChanIDs.FrequencyInHertz))...
                    'Hz_THR' num2str(tmpChanIDs.LowerLfpThreshold) '-'...
                    num2str(tmpChanIDs.UpperLfpThreshold) '_AVG'...
                    num2str(round(tmpChanIDs.AveragingDurationInMilliSeconds)) 'ms'];
                stimchannels{1,rl} = [stimName{rl} num2str(tmpChanIDs.RateInHertz) 'Hz_'...
                    num2str(tmpChanIDs.PulseWidthInMicroSecond) 'us'];
            else
                lfpsettings{rl,1}='LFP n/a';
                stimchannels{1,rl} = 'STIM n/a';
            end
            
        end
        
        
        
    case 2
        
        
end




end

