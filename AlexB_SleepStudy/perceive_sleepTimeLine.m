function [outMAT] = perceive_sleepTimeLine(inPS)
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
% 1. COMPLEX NAT SCHEME
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
    case 1 % Timeline
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

        lfpDayProc = replace(extractAfter(extractBefore(lfpDayNamesPRE,'T'),'x'),'_','/');
        lfpDayPro2 = datetime(lfpDayProc,'InputFormat','yyyy/MM/dd');
        lfpDayNames = lfpDayNamesPRE(ismember(lfpDayPro2,actDATES2));

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
%         monCOL = reshape(monthS,numel(monthS),1);
%         dayCOL = reshape(dayS,numel(dayS),1);
        hourCOL = reshape(hourS,numel(hourS),1);
        minCOL = reshape(minuteS,numel(minuteS),1);
        LFPaCOL = reshape(LFPallfSm,numel(LFPallfSm),1);
%         stimCOL = reshape(stimAll,numel(stimAll),1);

        % Remove smooth added edge data
        nanMin = isnan(minCOL);
        LFPaCOL(nanMin) = 0;

        % REMOVE NATs %%%%%%%%%%%%%%%%%%%%%%%%%% DEAL WITH NAT
%         [lfpDateInter] = natAssess(actDAYcol);
        % FILL 2 144 %%%%%%%%%%%%%%%%%%%%%%%%%% DEAL WITH NAT
        if any(isnat(actDAYcol))
            [lfpDateInter] = fill144all(actDAYcol);
        else
            lfpDateInter =  actDAYcol;
        end

        allLFPdat = [lfpDateInter , monCOL , dayCOL , hourCOL , minCOL ,...
                     LFPaCOL , stimCOL];

        % Save lfpDateInter
        saveStage1 = [saveLOC , filesep , 'STAGE1'];
        if ~exist(saveStage1,'dir')
            mkdir(saveStage1)
        end
        cd(saveStage1)
        svSt1n = ['LFP_Timeline_',num2str(subID),'.mat'];
        save(svSt1n,'allLFPdat');

    case 2

        % Process Actigraphy
        % Convert into 10 minute bins
        % Loop through date / hour / 10 min bins

        lfpDateREF = lfpDateInter;

        nonNatInd = ~isnat(lfpDateREF);
        tmpDLstr = cell(size(lfpDateREF));
        
        tmpDLstr(nonNatInd) = cellstr(datestr(lfpDateREF(nonNatInd)));
        lfpDateREF = tmpDLstr;
        lfpDateREF(nonNatInd) = extractBefore(tmpDLstr(nonNatInd),' ');

        isEMv = cellfun(@(x) ~isempty(x), lfpDateREF, 'UniformOutput',true);
        uniMDate = unique(lfpDateREF(isEMv));
        lfpDateREF(~isEMv) = {''};

        activityMean = zeros(144,length(uniMDate));
        activitySTD = zeros(144,length(uniMDate));
        WLMean = zeros(144,length(uniMDate));
        WLSTD = zeros(144,length(uniMDate));
        RLMean = zeros(144,length(uniMDate));
        RLSTD = zeros(144,length(uniMDate));
        GLMean = zeros(144,length(uniMDate));
        GLSTD = zeros(144,length(uniMDate));
        BLMean = zeros(144,length(uniMDate));
        BLSTD = zeros(144,length(uniMDate));
        WakeFrac = zeros(144,length(uniMDate));
        MaxInter = cell(144,length(uniMDate));

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
            % ADD 12 hours to ACT to get into 24 not 12 for LFP
            hourCOLnc = num2cell(hourCOL);
            minCOLtx = cellfun(@(x) num2str(x) , minCOLnc,'UniformOutput',false);
            hourCOLtx  = cellfun(@(x) num2str(x) , hourCOLnc,'UniformOutput',false);

            for mmi = 1:length(minCOLtx)
                if length(minCOLtx{mmi}) == 1
                    minCOLtx{mmi} = ['0',minCOLtx{mmi}];
                end
            end

            hourTX12hr = hourCOLtx;
            %%%%% CODE TO CONVERT LFP from 24 to 12

            hourMINlfpnn = cellfun(@(x,y) [x, ':', y],...
                hourTX12hr(uniDATEind), minCOLtx(uniDATEind),'UniformOutput',false);

            % EXPERIMENTAL !!!!!!!!!!!!!@@@@@@%%%%%
            if length(hourMINlfpnn) < 144
                addLEN = 144 - length(hourMINlfpnn);
                hourMINlfp2 = [hourMINlfpnn ; repmat({'NaN:NaN'},addLEN,1)];
            elseif length(hourMINlfpnn) > 144
                hourMINlfp2 = unique(hourMINlfpnn);
                if length(hourMINlfp2) ~= 144
                    hourMINlfp2 = hourMINlfp2(~ismember(hourMINlfp2,'0:00'));
                end
            end

            % CONVERT ACT TIME from 12 to 24
            [actHOUR24h] = convertBT1224('ACT',actDATEtab.Time);

            % Create searchable time array for actigraphy
            lastLOC = cellfun(@(x) find(x == ':',1,'last'), actHOUR24h,'UniformOutput',true);
            actINDtime = extractBefore(actHOUR24h,lastLOC);

            % 1: date, 2: Time, 3: ActivityMean, 4: ActivitySTD, 5: WhitelightMean
            % 6: WhitelightSTD, 7: RedLightMean, 8: RedLightSTD, 9:
            % GreenlightMean, 10: GreenlightSTD, 11: BluelightMean, 12:
            % BlulightSTD, 13: SleepWakefrac, 14: IntervalstatusMax
            act10minBin = cell(length(hourMINlfp2),14);
            act10minBinAll = cell(length(hourMINlfp2),1);

            % Loop through LFP times
            for lfpT = 1:length(hourMINlfp2)
                tmpLFPbin = hourMINlfp2{lfpT};
                % Remove blanks
                tmpLFPbin = replace(tmpLFPbin,' ','');
                % Start time bin
                binStart = find(matches(actINDtime,tmpLFPbin),1,'first');

                % If bin is not present because actigraphy was not recorded
                % I.E., off-wrist
                if isempty(binStart)
                    act10minBin{lfpT,3} = NaN;
                    act10minBin{lfpT,4} = NaN;
                    act10minBin{lfpT,5} = NaN;
                    act10minBin{lfpT,6} = NaN;
                    act10minBin{lfpT,7} = NaN;
                    act10minBin{lfpT,8} = NaN;
                    act10minBin{lfpT,9} = NaN;
                    act10minBin{lfpT,10} = NaN;
                    act10minBin{lfpT,11} = NaN;
                    act10minBin{lfpT,12} = NaN;
                    act10minBin{lfpT,13} = NaN;
                    continue
                else

                    binSIZE = (2*10) - 1; % Bins per min * num of minutes
                    binEnd = binStart + binSIZE;

                    if binEnd > length(actINDtime)
                        binEnd = length(actINDtime);
                    end

                    binTable = actDATEtab(binStart:binEnd,:);
                    act10minBinAll{lfpT} = binTable;

                    % date
                    act10minBin{lfpT,1} = tmpDI;
                    % time for 10 minute bin
                    act10minBin{lfpT,2} = actHOUR24h{binStart};
                    % average activity
                    tmpC1 = cellfun(@(x) str2double(x),...
                        binTable.Activity , 'UniformOutput',true);
                    act10minBin{lfpT,3} = mean(tmpC1(~isnan(tmpC1)));
                    act10minBin{lfpT,4} = std(tmpC1(~isnan(tmpC1)));
                    % whitelightM
                    tmpC2 = cellfun(@(x) str2double(x),...
                        binTable.("White Light") , 'UniformOutput',true);
                    act10minBin{lfpT,5} = mean(tmpC2(~isnan(tmpC2)));
                    act10minBin{lfpT,6} = std(tmpC2(~isnan(tmpC2)));
                    % redlightM
                    tmpC3 = cellfun(@(x) str2double(x),...
                        binTable.("Red Light") , 'UniformOutput',true);
                    act10minBin{lfpT,7} = mean(tmpC3(~isnan(tmpC3)));
                    act10minBin{lfpT,8} = std(tmpC3(~isnan(tmpC3)));
                    % greenlightM
                    tmpC4 = cellfun(@(x) str2double(x),...
                        binTable.("Green Light") , 'UniformOutput',true);
                    act10minBin{lfpT,9} = mean(tmpC4(~isnan(tmpC4)));
                    act10minBin{lfpT,10} = std(tmpC4(~isnan(tmpC4)));
                    % bluelightM
                    tmpC5 = cellfun(@(x) str2double(x),...
                        binTable.("Blue Light") , 'UniformOutput',true);
                    act10minBin{lfpT,11} = mean(tmpC5(~isnan(tmpC5)));
                    act10minBin{lfpT,12} = std(tmpC5(~isnan(tmpC5)));
                    % WakeFrac
                    sleepwakeT = cellfun(@(x) str2double(x),...
                        binTable.("Sleep/Wake") , 'UniformOutput',true);
                    act10minBin{lfpT,13} = sum(sleepwakeT)/length(sleepwakeT);
                    % IntervalstatusMax
                    uniStates = unique(binTable.("Interval Status"));
                    numStates = zeros(length(uniStates));
                    for ni = 1:length(uniStates)
                        numStates(ni) = sum(matches(binTable.("Interval Status"),uniStates{ni}));
                    end
                    [~,maxI] = max(numStates);
                    maxSTATE = uniStates{maxI};

                    act10minBin{lfpT,14} = maxSTATE;

                end % End of if statement checking bin presence
            end % End of bin loop
            % Finalize table
            activityMean(:,di) = cell2mat(act10minBin(:,3));
            activitySTD(:,di) = cell2mat(act10minBin(:,4));
            WLMean(:,di) = cell2mat(act10minBin(:,5));
            WLSTD(:,di) = cell2mat(act10minBin(:,6));
            RLMean(:,di) = cell2mat(act10minBin(:,7));
            RLSTD(:,di) = cell2mat(act10minBin(:,8));
            GLMean(:,di) = cell2mat(act10minBin(:,9));
            GLSTD(:,di) = cell2mat(act10minBin(:,10));
            BLMean(:,di) = cell2mat(act10minBin(:,11));
            BLSTD(:,di) = cell2mat(act10minBin(:,12));
            WakeFrac(:,di) = cell2mat(act10minBin(:,13));
            MaxInter(:,di) = act10minBin(:,14);

        end % End of date loop

        cd(saveLOC)

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
        outMAT.ActMean = activityMean;
        outMAT.ActSTD = activitySTD;
        outMAT.WLMean = WLMean;
        outMAT.WLSTD = WLSTD;
        outMAT.RLMean = RLMean;
        outMAT.RLSTD = RLSTD;
        outMAT.GLMean = GLMean;
        outMAT.GLSTD = GLSTD;
        outMAT.BLMean = BLMean;
        outMAT.BLSTD = BLSTD;
        outMAT.WakeFrac = WakeFrac;
        outMAT.MaxState = MaxInter;

        if matches(inPS.hemiS,'L')
            fileNAMEm = ['SPPD',num2str(inPS.subID),'_L_TimeLine.mat'];
        else
            fileNAMEm = ['SPPD',num2str(inPS.subID),'_R_TimeLine.mat'];
        end
        save(fileNAMEm,'outMAT');


end



end % End of Function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [conHourVec] = convertBT1224(time2con, convertVec)

if matches(time2con,'LFP')
    % converting from 24 to 12
    hourTX2nm = cellfun(@(x) str2double(x), convertVec, 'UniformOutput', true);
    hourTX2nm(hourTX2nm < 1) = abs(hourTX2nm(hourTX2nm < 1) - 12);
    hourTX2nm(hourTX2nm > 12) = abs(hourTX2nm(hourTX2nm > 12) - 12);
    conHourVec = cellstr(num2str(hourTX2nm));

elseif matches(time2con,'ACT')
    % converting from 12 to 24
    pm = contains(convertVec,'PM');
    am = contains(convertVec,'AM');
    firstHourloc = cellfun(@(x) find(x == ':',1,'first'), convertVec,'UniformOutput',true);
    getHOUR = extractBefore(convertVec,firstHourloc);
    hourNUM = cellfun(@(x) str2double(x) , getHOUR ,'UniformOutput',true);
    hourNUM(pm & hourNUM < 12) = hourNUM(pm & hourNUM < 12) + 12;
    hourNUM(am & hourNUM == 12) = 0;
    hourCeLL = cellfun(@(x) num2str(x) , num2cell(hourNUM), 'UniformOutput',false);
    indexCeLL = num2cell(firstHourloc);
    minSEcs = cellfun(@(x,y) extractAfter(x,y), convertVec, indexCeLL, 'UniformOutput', false);
    conHourVec = cellfun(@(x,z) [z,':',x], minSEcs,  hourCeLL , 'UniformOutput', false);
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [lfpDateInterALL] = natAssess(inACT2LFPbins)

% FIRST TRY KEEPING intake less than 144

% NEW PLAN FILL ALL TO 144 - FILL GAPs with NAT

% uniBINS = unique(inACT2LFPbins);
% dates2loop = uniqueDATES;

allDATES = day(inACT2LFPbins(~isnat(inACT2LFPbins)));
uniDATES = unique(day(inACT2LFPbins(~isnat(inACT2LFPbins))));
numDates = length(uniDATES);

lfpDateInterALL = NaT(size(inACT2LFPbins));
for ti = 1:numDates

    tmpDl = inACT2LFPbins(allDATES == uniDATES(ti));
    natLocs = find(isnat(tmpDl));
    natDiff = diff(find(isnat(tmpDl)));

    if isempty(natLocs) % if there are no NATS
        tmpDl3 = tmpDl;
        lfpDateInter = tmpDl3;
    elseif natLocs(1) == 1 % If first element is NAT
        firstNATs = find(natDiff ~= 1,1,'first');
        tmpDl2 = tmpDl(firstNATs+1:end);

        natLOCSm = find(isnat(tmpDl2));
        natDiffm = diff(natLOCSm);

        toRem = zeros(length(natLOCSm),1,'logical');
        for nii = 1:length(natDiffm)
            if natDiffm(nii) == 1
                toRem([nii,(nii + 1)]) = 1;
            else
                continue
            end
        end

        remLOG = ones(size(tmpDl2),"logical");
        remLOG(natLOCSm(toRem)) = false;
        tmpDl3 = tmpDl2(remLOG);

        % Interpolate single missing points
        natLogic = isnat(tmpDl3);
        natInd = find(natLogic);
        lfpDateInter = tmpDl3;
        for lli = 1:length(natInd)
            tmpIND = natInd(lli);
            preIND = tmpIND - 1;
            lfpDateInter(tmpIND) = lfpDateInter(preIND) + minutes(10);
        end

    elseif ~isempty(find(isnat(tmpDl), 1))  % If first element is not NAT and there are NATs
        natLOCSm = find(isnat(tmpDl));
        natDiffm = diff(natLOCSm);

        toRem = zeros(length(natLOCSm),1,'logical');
        for nii = 1:length(natDiffm)
            if natDiffm(nii) == 1
                toRem([nii,(nii + 1)]) = 1;
            else
                continue
            end
        end

        remLOG = ones(size(tmpDl),"logical");
        remLOG(natLOCSm(toRem)) = false;
        tmpDl3 = tmpDl(remLOG);

        % Interpolate single missing points
        natLogic = isnat(tmpDl3);
        if ~any(natLogic)
            lfpDateInter = tmpDl3;
        else
            natInd = find(natLogic);
            lfpDateInter = tmpDl3;
            for lli = 1:length(natInd)
                tmpIND = natInd(lli);
                preIND = tmpIND - 1;
                lfpDateInter(tmpIND) = lfpDateInter(preIND) + minutes(10);
            end
        end
    end
    % Get Length of lfpDateInter
    % Get first Nat of ALL
    % Fill
    tmpLen = length(lfpDateInter);
    natFirst = find(isnat(lfpDateInterALL),1,'first');
    lfpDateInterALL(natFirst:natFirst+tmpLen-1) = lfpDateInter;
end

lfpDateInterALL = lfpDateInterALL(~isnat(lfpDateInterALL));

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [out144fil] = fill144all(injagged)

allDATES = day(injagged(~isnat(injagged)));
uniDATES = unique(day(injagged(~isnat(injagged))));
numDates = length(uniDATES);
start = 1;
stop = 144;
out144fil = NaT(numDates*144,1);
for ti = 1:numDates
    tmpDl = injagged(allDATES == uniDATES(ti));
    nonnat = tmpDl(~isnat(tmpDl));
    nTest = nonnat(1);
    basicTimV = transpose(nTest:minutes(10):nTest + (hours(24) - minutes(10)));
    basicTimS = cellstr(datestr(basicTimV));

    tmpEx = extractAfter(basicTimS,' ');
    lastEl = cellfun(@(x) find(x == ':',1',"last"), tmpEx,...
        'UniformOutput', true);
    tmpEx2 = extractBefore(tmpEx,lastEl - 1);

    % Clean up tmpDl
    nonNatInd = ~isnat(tmpDl);
    tmpDLstr = cell(size(tmpDl));
    tmpDLstr(nonNatInd) = cellstr(datestr(tmpDl(nonNatInd)));
    tmpExDL = extractAfter(tmpDLstr(nonNatInd),' ');
    lastdlE = cellfun(@(x) find(x == ':',1',"last"), tmpExDL,...
        'UniformOutput', true);
    tmpdlEx2 = extractBefore(tmpExDL,lastdlE - 1);

    % Create file
    basicTimV(~ismember(tmpEx2,tmpdlEx2)) = NaT;

    out144fil(start:stop) = basicTimV;
    start = stop + 1;
    stop = (144 - 1) + start;
    
    

end






end