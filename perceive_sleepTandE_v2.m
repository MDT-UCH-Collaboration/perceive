function perceive_sleepTandE_v2(inPS)
% Github https://github.com/MDT-UCH-Collaboration

arguments
    inPS.userPC (1,1) string = "JATwork"
    inPS.subID (1,1) double = 1 % IN USE for CASE NUMBER
    inPS.postN (1,1) double = 1
    inPS.userDIR (1,1) string = "NA"
    inPS.seluDIR (1,1) logical = true
    inPS.saveDIR (1,1) string = "NA"
    inPS.selsDIR (1,1) logical = true
    inPS.stagE (1,1) double = 1
    inPS.studY (1,:) char = '20-2508'
    inPS.pltCl (1,1) logical = 0
end

%% OUTPUT
% The script generates BIDS inspired subject and session folders with the

%% TODO:
% 1. Remove less than zero from final LFP output
% 2. Figure what to output in CSV - long view (vs wide view)
% 

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
sessionFields = {'SessionDate','SessionEndDate','PatientInformation'};

cd(fileDIR)
initDir = dir('*.json');
jsonFiles = {initDir.name};

[indx,~] = listdlg('PromptString',{'Select a file.',...
    'Only one file can be selected at a time.',''},...
    'SelectionMode','single','ListString',jsonFiles);

json2load = jsonFiles{indx};

js = jsondecode(fileread(json2load));

% Tzoffset
tzOFF = js.ProgrammerUtcOffset;
tzOffs = strsplit(tzOFF,':');
tzOffsN = str2double(tzOffs{1});

switch inPS.stagE
    case 1 % Events
        infoFields = {'PatientEvents','EventSummary'};
        
        
    case 2 % Timeline
        infoFields = {'DiagnosticData'};
        
        dataOfInterest = js.(infoFields{1});
        
        lfpDAys = dataOfInterest.LFPTrendLogs.HemisphereLocationDef_Left;
        
        lfpDayNames = fieldnames(lfpDAys);
        
        monthS = zeros(144,length(lfpDayNames));
        dayS = zeros(144,length(lfpDayNames));
        hourS = nan(144,length(lfpDayNames));
        minuteS = nan(144,length(lfpDayNames));
        actDAYtm = NaT(144,length(lfpDayNames));
        LFPall = zeros(144,length(lfpDayNames));
        stimAll = zeros(144,length(lfpDayNames));
        
        for li = 1:length(lfpDayNames)
            
            tLFP = transpose(fliplr([lfpDAys.(lfpDayNames{li}).LFP]));
            tstim_mA = transpose(fliplr([lfpDAys.(lfpDayNames{li}).AmplitudeInMilliAmps]));
            timeD = transpose(fliplr({lfpDAys.(lfpDayNames{li}).DateTime}));
            
            [monthOI,dayOI,hourOI,minuteOI,actDayT] = getDT(timeD , tzOffsN);

            % convert minute column to floor round
            minuteOIc = floor(minuteOI/10)*10;
            % combine hour , minute , second
            durFind = duration(hourOI,minuteOIc,zeros(length(minuteOIc),1));
            
            % search for where to align times;
            [alignIND] = alignTime(durFind);
            
            monthS(alignIND,li) = monthOI;
            dayS(alignIND,li) = dayOI;
            hourS(alignIND,li) = hourOI;
            minuteS(alignIND,li) = minuteOI;
            actDAYtm(alignIND,li) = actDayT;
            LFPall(alignIND,li) = tLFP;
            stimAll(alignIND,li) = tstim_mA;

        end

        % Save out CSV file with timeline data
        % Month, Day, Hour, Minute, LFP mag, stimMA, actDayTime
        LFPallf = LFPall;
        LFPallf(LFPall > 7000) = 7000;
        
        % Fix outliers by removing and smooth
        [LFPallfm] = fixOutSm(LFPallf);
        
        % Smooth
        LFPallfSm = smoothdata(LFPallfm,1,'sgolay',6);
        
        % OutTable
        outTable = table(actDAYtm,monthS,dayS,hourS,minuteS,LFPallfSm,stimAll,...
            'VariableNames',{'FullDate','Month','Day','Hour','Minute',...
            'LFP_Mag','Stim_mA'});
        
        
        cd(saveLOC)
        fileNAME = ['SPPD',num2str(inPS.subID),'_TimeLine.csv'];
        writetable(outTable,fileNAME);

        
end



end % End of Function




function [monthOI,dayOI,hourOI,minuteOI,actDayT] = getDT(timeVEC , tzOFFtm)

monthOI = zeros(length(timeVEC),1);
dayOI = zeros(length(timeVEC),1);
hourOI = zeros(length(timeVEC),1);
minuteOI = zeros(length(timeVEC),1);
actDayT = NaT(length(timeVEC),1);

for ti = 1:length(timeVEC)
    
    tmpT = timeVEC{ti};
    tmpTT = datetime(replace(tmpT,{'T','Z'},{' ',''}));
    
    actDayT(ti) = tmpTT + hours(tzOFFtm);
    monthOI(ti) = actDayT(ti).Month;
    dayOI(ti) = actDayT(ti).Day;
    hourOI(ti) = actDayT(ti).Hour;
    minuteOI(ti) = actDayT(ti).Minute;
    
end

end



function [alignIND] = alignTime(inTIME)

AMblock1 = duration(6,00,0);
AMblock2 = duration(23,50,0);

PMblock1 = duration(0,00,0);
PMblock2 = duration(5,50,0);

amBlock = linspace(AMblock1,AMblock2,108);
pmBlock = linspace(PMblock1,PMblock2,36);

% Start at 6 AM
allBlock = [transpose(amBlock) ; transpose(pmBlock)];

% Search through inTIME and line up with allBlock indicies in alignIND
% Check on the 143 and 50 file.

alignIND = zeros(size(inTIME));
for iT = 1:length(inTIME)
    
    % Input time row
    tTime = inTIME(iT);
    
    % Find where located in allBlock
    tIND  = find(ismember(allBlock,tTime));
    
    % Store in alignIND
    alignIND(iT) = tIND;
    
end

end


function [outSmooth] = fixOutSm(inUNvec)


outSmooth = inUNvec;
for si = 1:size(inUNvec,2)
    
    tmpVEC = inUNvec(:,si);
    tmpVECz = tmpVEC(tmpVEC ~= 0);
    
    tmpM = mean(tmpVECz);
    tmpS = std(tmpVECz);
    
    tmpTH = tmpM + (tmpS*2);
    
    abovEInd = find(tmpVEC > tmpTH);
    allNUMS = 1:144;
    remAbov = ~ismember(allNUMS,abovEInd);

    for ai = 1:length(abovEInd)
        tmpAi = abovEInd(ai);
        bef = tmpAi - 6:tmpAi - 1;
        befD = bef(remAbov(bef));
        
        aft = tmpAi + 1:tmpAi + 6;
        aftD = aft(remAbov(aft));
        
        menBlk = mean(tmpVEC([befD , aftD]));
        tmpVEC(tmpAi) = menBlk;
    end
    outSmooth(:,si) = tmpVEC;
  
end

end

