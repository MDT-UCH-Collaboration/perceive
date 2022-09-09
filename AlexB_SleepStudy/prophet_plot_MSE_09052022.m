cd('E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\prophetOUT')

% rawData = readtable("SPPD3_L_Prophet.csv");
% forDdata = readtable("SPPD3_L_forcast.csv");

csvDir = dir('*.csv');
csvDirn = {csvDir.name};

csvParts = split(csvDirn,'_');
sppdComb = cellfun(@(x,y) [x,'_',y],csvParts(:,:,1),csvParts(:,:,2),'UniformOutput',false );
conditionS = csvParts(:,:,3);

sppdUnique = unique(sppdComb);

for si = 1:length(sppdUnique)

    sppdIND = matches(sppdComb, sppdUnique{si});

    sppINDc = conditionS(sppdIND);

    for ci = 1:2
        tmpFC = sppINDc{ci};
        if contains(tmpFC,'Prophet')
            rawFileN = [sppdUnique{si},'_',sppINDc{ci}];
            rawData = readtable(rawFileN);
        else
            forFileN = [sppdUnique{si},'_',sppINDc{ci}];
            forDdata = readtable(forFileN);
        end

    end

    % MSE
    sharedINdsFOR = ismember(forDdata.ds,rawData.ds);
    actual = rawData.y;
    prediCTEd = forDdata.yhat(sharedINdsFOR);
    mseOUT = mse(actual, prediCTEd);


    % Plot
    figure;
    plot(rawData.ds,rawData.y,'k.','MarkerSize',15)

    hold on

    plot(forDdata.ds,forDdata.yhat,'r-','MarkerSize',5)

    yPAT = [transpose(forDdata.yhat_lower) , fliplr(transpose(forDdata.yhat_upper))];
    xPAT = [transpose(forDdata.ds) , fliplr(transpose(forDdata.ds))];

    p = patch(xPAT,yPAT,'r');
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.2;

end


%%

% plot(rawData.ds,rawData.y,'k.','MarkerSize',15)
%
% hold on
%
% plot(forDdata.ds,forDdata.yhat,'r-','MarkerSize',5)
%
% yPAT = [transpose(forDdata.yhat_lower) , fliplr(transpose(forDdata.yhat_upper))];
% xPAT = [transpose(forDdata.ds) , fliplr(transpose(forDdata.ds))];
%
% p = patch(xPAT,yPAT,'r');
% p.EdgeColor = 'none';
% p.FaceAlpha = 0.2;


cd('E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\prophetOUT')

% rawData = readtable("SPPD3_L_Prophet.csv");
% forDdata = readtable("SPPD3_L_forcast.csv");

csvDir = dir('*.csv');
csvDirn = {csvDir.name};

csvParts = split(csvDirn,'_');
sppdComb = cellfun(@(x,y) [x,'_',y],csvParts(:,:,1),csvParts(:,:,2),'UniformOutput',false );
conditionS = csvParts(:,:,3);

sppdUnique = unique(sppdComb);

mseOUT = zeros(1,length(sppdUnique));
for si = 1:length(sppdUnique)

    sppdIND = matches(sppdComb, sppdUnique{si});

    sppINDc = conditionS(sppdIND);

    for ci = 1:2
        tmpFC = sppINDc{ci};
        if contains(tmpFC,'Prophet')
            rawFileN = [sppdUnique{si},'_',sppINDc{ci}];
            rawData = readtable(rawFileN);
        else
            forFileN = [sppdUnique{si},'_',sppINDc{ci}];
            forDdata = readtable(forFileN);
        end

    end

    % MSE
    sharedINdsFOR = ismember(forDdata.ds,rawData.ds);
    actual = rawData.y;
    prediCTEd = forDdata.yhat(sharedINdsFOR);
    mseOUT(si) = mse(actual, prediCTEd);

end


