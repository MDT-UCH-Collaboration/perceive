function [outMAT] = perceive_sleepTandE_v3(inPS)
% Github https://github.com/MDT-UCH-Collaboration

arguments
    inPS.userPC (1,1) string = "JATwork"
    inPS.subID (1,1) double = 1 % IN USE for CASE NUMBER
    inPS.postN (1,1) double = 1
    inPS.userDIR (1,1) string = "NA"
    inPS.hemiS (1,1) string = "L"
    inPS.seluDIR (1,1) logical = true
    inPS.saveDIR (1,1) string = "NA"
    inPS.actDIR (1,1) string = "NA"
    inPS.selsDIR (1,1) logical = true
    inPS.stagE (1,1) double = 1
    inPS.studY (1,:) char = '20-2508'
    inPS.pltCl (1,1) logical = 0
    inPS.tabLOC (1,1) string = "NA"
    inPS.overSAT (1,1) logical = 1
    inPS.jsonDAT (1,1) string = "NA"
end

%% OUTPUT
% 

%% TODO:
% 1. Fix Event time
%

if inPS.seluDIR && strcmp(inPS.userDIR,"NA")
    [fileDIR] = uigetdir();
else
    if matches(inPS.hemiS,'L')
        fileDIR = [char(inPS.userDIR) , filesep , 'Left'];
    else
        fileDIR = [char(inPS.userDIR) , filesep , 'Right'];
    end
end

if strcmp(inPS.actDIR, "NA")
    [actLOC] = uigetdir();
else
    actLOC = [char(inPS.actDIR) , num2str(inPS.subID),...
        '\ACT_data\Summary'];
end

if inPS.selsDIR && strcmp(inPS.saveDIR,"NA")
    [saveLOC] = uigetdir();
else
    saveLOC = inPS.saveDIR;
end

if strcmp(inPS.tabLOC,"NA")
    [tabFname] = uigetfile();
    dataTABLE = readtable(tabFname);
else
    dataTABLE = readtable(inPS.tabLOC);
end

patTable = dataTABLE(inPS.subID,:);

sessionFields = {'SessionDate','SessionEndDate','PatientInformation'};

cd(fileDIR)

if strcmp(inPS.jsonDAT,"NA")

    initDir = dir('*.json');
    jsonFiles = {initDir.name};

    [indx,~] = listdlg('PromptString',{'Select a file.',...
        'Only one file can be selected at a time.',''},...
        'SelectionMode','single','ListString',jsonFiles);

    json2load = jsonFiles{indx};
else
    json2load = inPS.jsonDAT;
end

js = jsondecode(fileread(json2load));

% Tzoffset
tzOFF = js.ProgrammerUtcOffset;
tzOffs = strsplit(tzOFF,':');
tzOffsN = str2double(tzOffs{1});

switch inPS.stagE
    case 1 % Events
        infoFields = {'PatientEvents','EventSummary','DiagnosticData'};
        js.DiagnosticData.LfpFrequencySnapshotEvents
        dataOfInterest = js.(infoFields{3});
        
        % Struct
        % Level 1
        % - Type of event
        % Level 2
        % -- Event start time
        % -- Matrix of FFT [column = event #, row = sample #]
        % -- Matrix of STim [column = event #, row = sample #]
        % -- Matrix of LFPmag [column = event #, row = sample #]
        
        allEvents = cell(length(dataOfInterest.LfpFrequencySnapshotEvents),1);
        
        for ae = 1:length(allEvents)
            allEvents{ae} = dataOfInterest.LfpFrequencySnapshotEvents{ae,1}.EventName;
        end
        
        allEventsS = cellfun(@(x) replace(x,' ',''), allEvents, 'UniformOutput',false);
        
        uniEVENTS = unique(allEventsS);
        
        outEVENTS = struct;
        for oi = 1:length(uniEVENTS)
            outEVENTS.(uniEVENTS{oi}).count = 1;
        end
        
        for di = 1:length(dataOfInterest.LfpFrequencySnapshotEvents)
            
            tmpDat = dataOfInterest.LfpFrequencySnapshotEvents{di,1};
            
            tEvent = tmpDat.EventName;
            tEventS = replace(tEvent,' ','');
            tIND = ismember(uniEVENTS,tEventS);
            
            inCount = outEVENTS.(uniEVENTS{tIND}).count;
            
            if ~isfield(tmpDat,'LfpFrequencySnapshotEvents')
                continue
            end
            
            [monthOI,dayOI,hourOI,minuteOI,actTIMEDay] = getDT({tmpDat.DateTime} , tzOffsN);
            
            outEVENTS.(uniEVENTS{tIND}).DateTime{inCount,1} = actTIMEDay;

            outEVENTS.(uniEVENTS{tIND}).Month(inCount,1) = monthOI;
            outEVENTS.(uniEVENTS{tIND}).Day(inCount,1) = dayOI;
            outEVENTS.(uniEVENTS{tIND}).Hour(inCount,1) = hourOI;
            outEVENTS.(uniEVENTS{tIND}).Minute(inCount,1) = minuteOI;
            
            lfpDAT = tmpDat.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left;

            % Smooth data
            lfpDATz = lfpDAT.FFTBinData;
            lfpDATzSM = smoothdata(lfpDATz,'sgolay',5);
            % Remove 0s
            lfpDATzRM = lfpDATzSM;
            lfpDATzRM(lfpDATzSM < 0) = 0;
            
            outEVENTS.(uniEVENTS{tIND}).FFTBinData(:,inCount) = lfpDATzRM;
            
            outEVENTS.(uniEVENTS{tIND}).Frequency(:,inCount) = lfpDAT.Frequency;
            
            outEVENTS.(uniEVENTS{tIND}).count = inCount + 1;
            
        end
        
        % Clean up struct
        for ci = 1:length(uniEVENTS)
            outEVENTS.(uniEVENTS{ci}) = rmfield(outEVENTS.(uniEVENTS{ci}),'count');
        end
        
        % Save mat file
        cd(saveLOC)
        fileNAMEm = ['SPPD',num2str(inPS.subID),'_Events.mat'];
        save(fileNAMEm,'outEVENTS');
        
        % Save Excel file 
        outTableE = table;
        for ci = 1:length(uniEVENTS)
            % Time
            timeAll = cell(100,length(outEVENTS.(uniEVENTS{ci}).DateTime));
            monthOI = zeros(100,length(outEVENTS.(uniEVENTS{ci}).DateTime));
            dayOI = zeros(100,length(outEVENTS.(uniEVENTS{ci}).DateTime));
            hourOI = zeros(100,length(outEVENTS.(uniEVENTS{ci}).DateTime));
            minuteOI = zeros(100,length(outEVENTS.(uniEVENTS{ci}).DateTime));
            for ti = 1:length(outEVENTS.(uniEVENTS{ci}).DateTime)
               timeAll(1:100,ti) = repmat(outEVENTS.(uniEVENTS{ci}).DateTime(ti),100,1);
               monthOI(1:100,ti) = repmat(outEVENTS.(uniEVENTS{ci}).Month(ti),100,1); 
               dayOI(1:100,ti) = repmat(outEVENTS.(uniEVENTS{ci}).Day(ti),100,1);
               hourOI(1:100,ti) = repmat(outEVENTS.(uniEVENTS{ci}).Hour(ti),100,1);
               minuteOI(1:100,ti) = repmat(outEVENTS.(uniEVENTS{ci}).Minute(ti),100,1);
            end
            timeAll2 = timeAll(:);
            monthOI2 = monthOI(:);
            dayOI2 = dayOI(:);
            hourOI2 = hourOI(:);
            minuteOI2 = minuteOI(:);
            
            outTableE.([uniEVENTS{ci},'Time']) = timeAll2;
            outTableE.([uniEVENTS{ci},'_Month']) = monthOI2;
            outTableE.([uniEVENTS{ci},'_Day']) = dayOI2;
            outTableE.([uniEVENTS{ci},'_Hour']) = hourOI2;
            outTableE.([uniEVENTS{ci},'_Minute']) = minuteOI2;
            % FFT
            outTableE.([uniEVENTS{ci},'FFT']) = outEVENTS.(uniEVENTS{ci}).FFTBinData(:);
            % Frequency
            outTableE.([uniEVENTS{ci},'Freq']) = outEVENTS.(uniEVENTS{ci}).Frequency(:);
        end
        
        cd(saveLOC)
        fileNAME = ['SPPD',num2str(inPS.subID),'_Events.csv'];
        writetable(outTableE,fileNAME);
        
        
    case 2 % Timeline
        infoFields = {'DiagnosticData'};

        groupID = [js.Groups.Final.ActiveGroup];
        activeGROUP = js.Groups.Final(groupID).ProgramSettings.SensingChannel;
        activeSenseChan = activeGROUP.Channel;
        activeSenseFreq = activeGROUP.SensingSetup.FrequencyInHertz;
        
        dataOfInterest = js.(infoFields{1});

        cd(actLOC)
        load(['SPPD',num2str(inPS.subID) ,'_ACT_DATA.mat'],'dataTable');

        % First first row of non nan
%         allNANs = cellfun(@(x) ~matches(x,'NaN'), dataTable.("White Light"), 'UniformOutput', true);
        onWristData = cellfun(@(x) str2double(x), dataTable.("Off-Wrist Status"));
        actTABLE = dataTable(~onWristData,:);

        % Extract Unique dates
        actDATES = unique(actTABLE.Date);
        actDATES2 = sortrows(datetime(actDATES,'InputFormat','M/dd/yyyy'));

        if contains(inPS.hemiS,"L")
            lfpDAys = dataOfInterest.LFPTrendLogs.HemisphereLocationDef_Left;
        else
            lfpDAys = dataOfInterest.LFPTrendLogs.HemisphereLocationDef_Right;
        end
        
        lfpDayNamesPRE = fieldnames(lfpDAys);

        % CHECK FOR AND FIND START AND END DATES - WITH DIARY DATA
        %%%%% FUTURE - MAKE AN input argument to switch between diary and
        %%%%% actigraphy
%         dstart = patTable.diaryStart;
%         dstartF = extractBefore(string(dstart),6);
% 
%         dend = patTable.diaryEnd;
%         dendF = extractBefore(string(dend),6);

        lfpDayProc = replace(extractAfter(extractBefore(lfpDayNamesPRE,'T'),'x'),'_','/');
        lfpDayPro2 = datetime(lfpDayProc,'InputFormat','yyyy/MM/dd');

        actLFPdates = lfpDayPro2(ismember(lfpDayPro2,actDATES2));
        lfpDayNames = lfpDayNamesPRE(ismember(lfpDayPro2,actDATES2));

%         dstartIND = find(contains(lfpDayProc,dstartF));
%         dendIND = find(contains(lfpDayProc,dendF));

%         lfpDayNames = lfpDayNamesPRE(dstartIND:dendIND);
       
        monthS = zeros(144,length(lfpDayNames));
        dayS = zeros(144,length(lfpDayNames));
        hourS = nan(144,length(lfpDayNames));
        minuteS = nan(144,length(lfpDayNames));
        actDAYtm = NaT(144,length(lfpDayNames));
        LFPall = zeros(144,length(lfpDayNames));
        stimAll = zeros(144,length(lfpDayNames));
        timXax = [];
        
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
            [alignIND , allBlok] = alignTime(durFind);
            
            monthS(alignIND,li) = monthOI;
            dayS(alignIND,li) = dayOI;
            hourS(alignIND,li) = hourOI;
            minuteS(alignIND,li) = minuteOI;
            actDAYtm(alignIND,li) = actDayT;
            LFPall(alignIND,li) = tLFP;
            stimAll(alignIND,li) = tstim_mA;
            timXax = allBlok;
            
        end
        
        % Save out CSV file with timeline data
        % Month, Day, Hour, Minute, LFP mag, stimMA, actDayTime
        LFPallf = LFPall;
        if inPS.overSAT
            LFPallf(LFPall > 7000) = 7000;
        else
            LFPallf(LFPall > 80000) = 80000;
        end

        % Fix outliers by removing and smooth
        [LFPallfm] = fixOutSm(LFPallf);
        
        % Smooth
        LFPallfSm = smoothdata(LFPallfm,1,'sgolay',6);
        
        % Tall Table
        actDAYcol = reshape(actDAYtm,numel(actDAYtm),1);
        monCOL = reshape(monthS,numel(monthS),1);
        dayCOL = reshape(dayS,numel(dayS),1);
        hourCOL = reshape(hourS,numel(hourS),1);
        minCOL = reshape(minuteS,numel(minuteS),1);
        %         minMapCOL = reshape(actDAYtm,numel(actDAYtm),1);
        LFPaCOL = reshape(LFPallfSm,numel(LFPallfSm),1);
        stimCOL = reshape(stimAll,numel(stimAll),1);
        
        %         timeCOL = repmat(minuteOIc,144*size(hourCOL,2),1);
        %         timeMAT = repmat(minuteOIc,1,size(hourCOL,2));
        
        % Remove smooth added edge data
        nanMin = isnan(minCOL);
        LFPaCOL(nanMin) = 0;
        
        % OutTable
        outTable = table(actDAYcol,monCOL,dayCOL,hourCOL,minCOL,LFPaCOL,stimCOL,...
            'VariableNames',{'FullDate','Month','Day','Hour','Minute',...
            'LFP_Mag','Stim_mA'});
        
        %         nanMin = isnan(outTable.Minute);
        %         outTable.LFP_Mag(nanMin) = 0;

        % Process Actigraphy
        % Convert into 10 minute bins
        % Loop through date / hour / 10 min bins
        lfpDateREF  = datetime(extractBefore(cellstr(datestr(actDAYcol)),' '));
        uniMDate = unique(lfpDateREF);

        for di = 1:length(uniMDate)
            % Unique date
            tmpDI = uniMDate(di);

            % temp act table
            dconvert = cellfun(@(x) replace(x,'/','-'),actTABLE.Date,"UniformOutput",false);
            actDATEn = datetime(dconvert,'InputFormat','M-dd-uuuu');

            % Logical date index - LFP
            uniDATEind = ismember(lfpDateREF,tmpDI);

            % Logical date index - actigraphy
            actDATEind = ismember(actDATEn, tmpDI);
            actDATEtab = actTABLE(actDATEind,:);

            % Create HOUR:MIN index from LFP
            % H:MM num2cell(hourCOL(uniDATEind)) num2cell(minCOL(uniDATEind))
            minCOLnc = num2cell(minCOL);
            % SUBTRACT 12 hours from LFP to get into 12 not 24 for ACT
            hourCOLnc = num2cell(hourCOL);
            minCOLtx = cellfun(@(x) num2str(x) , minCOLnc,'UniformOutput',false);
            hourCOLtx  = cellfun(@(x) num2str(x) , hourCOLnc,'UniformOutput',false);
            
            for mmi = 1:length(minCOLtx) 
                if length(minCOLtx{mmi}) == 1
                    minCOLtx{mmi} = ['0',minCOLtx{mmi}];
                end

            end

            hourMINlfp2 = cellfun(@(x,y) [x, ':', y],...
                hourCOLtx(uniDATEind), minCOLtx(uniDATEind),'UniformOutput',false);
            



        end

        
        cd(saveLOC)
        if matches(inPS.hemiS,'L')
            fileNAME = ['SPPD',num2str(inPS.subID),'_L_TimeLine.csv'];
        else
            fileNAME = ['SPPD',num2str(inPS.subID),'_R_TimeLine.csv'];
        end
        writetable(outTable,fileNAME);
        
        LFPaMAT = reshape(LFPaCOL,144,size(hourS,2));
        
        outMAT.actTime = actDAYtm;
        outMAT.month = monthS;
        outMAT.day = dayS;
        outMAT.hour = hourS;
        outMAT.minu = minuteS;
        outMAT.LFP = LFPaMAT;
        outMAT.Stim = stimAll;
        outMAT.TimeX = cellstr(datestr(timXax));
        outMAT.senseChan = activeSenseChan;
        outMAT.senseFreq = activeSenseFreq; 
        
        if matches(inPS.hemiS,'L')
            fileNAMEm = ['SPPD',num2str(inPS.subID),'_L_TimeLine.mat'];
        else
            fileNAMEm = ['SPPD',num2str(inPS.subID),'_R_TimeLine.mat'];
        end
        save(fileNAMEm,'outMAT');
        
        
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



function [alignIND , allBlock] = alignTime(inTIME)

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
        if tmpAi >= 138
            tmpAi = 138;
        end
        if tmpAi < 7
            tmpAi = 7;
        end
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


