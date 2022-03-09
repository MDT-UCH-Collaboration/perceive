function [] = findEventActtime()
% Clean up Event data with Actigraphy dates

% Change directory
name = getenv('COMPUTERNAME');
if matches(name,'DESKTOP-95LU6PO') % Home PC
    mainDIR = 'D:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data';
else
    mainDIR = 'C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data';
end

% 1. Loop through cases

cd(mainDIR)
dirS1 = dir();
dirS2 = {dirS1.name};
dirS3 = dirS2(~ismember(dirS2,{'.','..'}));
allCases = dirS3;

for ai = 1:length(allCases)

    tmpDir = [mainDIR , filesep , allCases{ai} , filesep , 'MAT_data'];
    cd(tmpDir)

    allMat = dir('*.mat');
    nmMat = {allMat.name};
    tlMat = nmMat{contains(nmMat,'TimeLine')};
    evMat = nmMat{contains(nmMat,'Events')};

    load(tlMat,'outMAT')
    load(evMat,'outTableE')

    % Do something
    allActTimes = outMAT.actTime;
    % Need one index per column
    % Day Matrix
    allDays = day(allActTimes(1,:));
    % Hour Matrix
    allHours = hour(allActTimes);
    % Minute Matrix in 10s
    allMins = floor(minute(allActTimes)/10)*10;

    % Loop through events
    fNames = fieldnames(outTableE);
    % Event indices
    eventINdices = zeros(length(fNames),length(allDays));
    for eii = 1:length(fNames)
        evTimes = outTableE.(fNames{eii}).Time(:);

        uniDays = unique(outTableE.(fNames{eii}).Day);
        for ui = 1:length(uniDays)

            uiTmp = uniDays(ui);
            tmpTime = evTimes{find(outTableE.(fNames{eii}).Day == uiTmp, 1, 'first')};
            tmpTimeDay = day(tmpTime);
            tmpTimeHour = hour(tmpTime);
            tmpTimeMin = floor(minute(tmpTime)/10)*10;

            % Find day column
            dayLogical = ismember(allDays(1,:),tmpTimeDay);
            if ~any(dayLogical)
                continue
            end
            dayHCol = allHours(:,dayLogical);
            dayMCol = allMins(:,dayLogical);

            % Find hour block and min index
            eventINdices(eii,dayLogical) =...
                find(dayHCol == tmpTimeHour & dayMCol == tmpTimeMin);
        end
    end
end

