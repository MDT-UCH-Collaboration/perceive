function [outTable] = ocd_dbs_rmanova()

% cd('D:\Dropbox\OCD_DBS_JSON')
% cd('D:\Dropbox\OCD_DBS_JSON\MP_CaseReport_2023')
% cd('C:\Users\johna\Dropbox\OCD_DBS_JSON\MP_CaseReport_2023')
cd('E:\Dropbox\OCD_DBS_JSON\MP_CaseReport_2023')

trialNall = [1 1 2 2];
hemisAll = {'L','R','L','R'};

% L Post C = 6 B = 6 
trialAllCells = cell(1,4);
for tni = 1:length(trialNall)

    tempTrial = trialNall(tni);
    tempHemi = hemisAll{tni};
    tempPrePos = 'after';

    [tmpfname] = getJSONofInt(tempTrial , tempHemi , tempPrePos);

    jsonfile = jsondecode(fileread(tmpfname));
    [procMat , tmpHzTRx , senseELall2] = getProcessMat(jsonfile);
    contactByBand = nan(5,6);
    for conI = 1:height(procMat)
        tmpContact = procMat(conI,:);
        % Extract max peak uVp and Hz
        bandsti = [1 , 4, 8, 13, 31];
        bandsto = [3 , 7, 12, 30 50];

        for bandi = 1:length(bandsti)

            bandINDEX = tmpHzTRx >= bandsti(bandi) & tmpHzTRx <= bandsto(bandi);
            tmpROWb = tmpContact(bandINDEX);
            bandMean = mean(tmpROWb);
            contactByBand(bandi,conI) = bandMean;

        end

    end

    conBandTab = array2table(contactByBand,"RowNames",{'delta','theta','alpha','beta','gamma'},...
        'VariableNames',senseELall2);

    trialAllCells{tni} = conBandTab;

end


allHemi = {};
allPow = [];
allBand = [];
allCont = {};
for tni = 1:length(trialNall)

    tempHemi = hemisAll(tni);

    powVals = reshape(table2array(trialAllCells{tni}),...
        numel(table2array(trialAllCells{tni})),1);

    tempBand = repmat(trialAllCells{tni}.Row,6,1);

    tempCont = reshape(repmat(trialAllCells{tni}.Properties.VariableNames,5,1),...
        numel(repmat(trialAllCells{tni}.Properties.VariableNames,5,1)),1);

    tempHemi2 = repmat(tempHemi,size(tempCont));

    allHemi = [allHemi ; tempHemi2];
    allPow = [allPow ; powVals];
    allBand = [allBand ; tempBand];
    allCont = [allCont ; tempCont];

end

outTable = table(allHemi, allPow , allBand , allCont,'VariableNames',...
    {'Hemi','Power','Band','Contact'});


end














function [jsonFnameTab] = getJSONofInt(trialN , hemisphere , prePOST)

flist = dir('*.json');
flistjson = {flist.name};

if trialN == 1
    triRows = flistjson(contains(flistjson , '1sttrial'));
else
    triRows = flistjson(contains(flistjson , '2ndtrial'));
end

if matches(hemisphere,'L')
    hemiRows = triRows(contains(triRows, 'Left'));
else
    hemiRows = triRows(contains(triRows, 'Right'));
end

if matches(prePOST,'after')
    hemiRows2 = hemiRows(contains(hemiRows, 'after'));
else
    hemiRows2 = hemiRows(contains(hemiRows, 'before'));
end

% fileName , hemisphere, trialN
jsonFnameTab = hemiRows2{1};

% jsonFnameTab = table(hemiRows3,repmat({num2str(trialN)},height(hemiRows3),1),...
%     repmat({hemisphere},height(hemiRows3),1),...
%     'VariableNames',{'jsonFiles','TrialN','HemiSide'});


end





function [procMat , tmpHzTRx , senseELall2] = getProcessMat(jsonfile)


tmpLFPbrainSense = jsonfile.LFPMontage;
tmpLFPtab = struct2table(tmpLFPbrainSense);


smothLFPtrim = zeros(6,72);
senseELall = cell(6,2);
for ti = 1:height(tmpLFPtab)

    senseEle = tmpLFPtab.SensingElectrodes{ti};
    senseEleName = extractAfter(senseEle,'.');
    [senseELnums] = translateNums(senseEleName);

    tmpLFP = tmpLFPtab.LFPMagnitude{ti};
    tmpHz = tmpLFPtab.LFPFrequency{ti};
    tmpHZcut = tmpHz < 70;

    tmpLFPsm = smoothdata(tmpLFP,'gaussian',7);
    tmpLFPtrim = tmpLFPsm(tmpHZcut);

    % plot(tmpHz(tmpHZcut),tmpLFPsm(tmpHZcut))
    smothLFPtrim(ti,:) = tmpLFPtrim;
    senseELall(ti,:) = senseELnums;

end

tmpHzTRx = tmpHz(tmpHZcut);

% Normalize
procMat = smothLFPtrim;
% unfurl = reshape(smothLFPtrim,numel(smothLFPtrim),1);
% normLFP = normalize(unfurl,'range');
% procMat = reshape(normLFP,6,72);
% plot(transpose(normSlfp))
% senseELallS = senseELall(high2low,:);
senseELall2 = join(senseELall,'-');

end



function [getNUms] = translateNums(eleName)

splitNames = split(eleName,'_');

getNUms = cell(1,2);
numCount = 1;
for si = 1:length(splitNames)

    if ~matches(splitNames{si},{'ZERO','ONE','TWO','THREE'})
        continue
    else
        switch splitNames{si}
            case 'ZERO'
                getNUms{numCount} = '0';
            case 'ONE'
                getNUms{numCount} = '1';
            case 'TWO'
                getNUms{numCount} = '2';
            case 'THREE'
                getNUms{numCount} = '3';
        end
        numCount = numCount + 1;
    end

end

end


function [outTABLE] = normalDAT(inTABLE)


allDATA = [inTABLE.LFPdata{1} ; inTABLE.LFPdata{2}];

unfurl = reshape(allDATA,numel(allDATA),1);
normLFP = normalize(unfurl,'range');
procMat = reshape(normLFP,12,72);

outTABLE = inTABLE;
outTABLE.LFPdata{1} = procMat(1:6,:);
outTABLE.LFPdata{2} = procMat(7:12,:);


end