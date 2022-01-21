function [outMAT] = perceive_sleepEvents(inPS)
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

% if strcmp(inPS.tabLOC,"NA")
%     [tabFname] = uigetfile();
%     dataTABLE = readtable(tabFname);
% else
%     dataTABLE = readtable(inPS.tabLOC);
% end
%
% patTable = dataTABLE(inPS.subID,:);

% sessionFields = {'SessionDate','SessionEndDate','PatientInformation'};

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

        if ~isfield(js.DiagnosticData,'LfpFrequencySnapshotEvents')
            return
        else

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
                try
                    allEvents{ae} = dataOfInterest.LfpFrequencySnapshotEvents{ae,1}.EventName;
                catch
                    allEvents{ae} = dataOfInterest.LfpFrequencySnapshotEvents(1).EventName;
                end
            end

            allEventsS = cellfun(@(x) replace(x,' ',''), allEvents, 'UniformOutput',false);

            uniEVENTS = unique(allEventsS);

            outEVENTS = struct;
            for oi = 1:length(uniEVENTS)
                outEVENTS.(uniEVENTS{oi}).count = 1;
            end

            for di = 1:length(dataOfInterest.LfpFrequencySnapshotEvents)

                try
                    tmpDat = dataOfInterest.LfpFrequencySnapshotEvents{di,1};
                catch
                    tmpDat = dataOfInterest.LfpFrequencySnapshotEvents(1);
                end


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

                hemiField = fieldnames(tmpDat.LfpFrequencySnapshotEvents);

                if contains(hemiField,'Left')
                    lfpDAT = tmpDat.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left;
                    hemiS = 'L';
                else
                    lfpDAT = tmpDat.LfpFrequencySnapshotEvents.HemisphereLocationDef_Right;
                    hemiS = 'R';
                end

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
            outTableE = struct;
            for ci = 1:length(uniEVENTS)
                % Time
                % Check if event has fields
                if isfield(outEVENTS.(uniEVENTS{ci}),'Hour')

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

                    outTableE.([uniEVENTS{ci}]).Time = timeAll2;
                    outTableE.([uniEVENTS{ci}]).Month = monthOI2;
                    outTableE.([uniEVENTS{ci}]).Day = dayOI2;
                    outTableE.([uniEVENTS{ci}]).Hour = hourOI2;
                    outTableE.([uniEVENTS{ci}]).Minute = minuteOI2;
                    % FFT
                    outTableE.([uniEVENTS{ci}]).FFT = outEVENTS.(uniEVENTS{ci}).FFTBinData(:);
                    % Frequency
                    outTableE.([uniEVENTS{ci}]).Freq = outEVENTS.(uniEVENTS{ci}).Frequency(:);
                else
                    continue
                end
            end

            cd(saveLOC)
            if matches(hemiS,'L')
                fileNAME = ['SPPD',num2str(inPS.subID),'_L_Events.mat'];
            else
                fileNAME = ['SPPD',num2str(inPS.subID),'_R_Events.mat'];
            end
            save(fileNAME,'outTableE');
        end
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
