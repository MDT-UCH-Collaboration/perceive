cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\SPPD2\MAT_data')


%%

load("SPPD2_L_TimeLine.mat")

%% 

daySu = unique(dataTable.Date);

allSleept = cell(length(daySu),1);
for di = 1:length(daySu)

    dayI = daySu{di};
    dayTable = dataTable(ismember(dataTable.Date,dayI),:);
    allTimes = dayTable.Time;
    sleepTimes = allTimes(contains(dayTable.("Interval Status"),'REST'));

    allSleept{di} = sleepTimes;

end

%%

