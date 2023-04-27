function [] = ocd_dbs_cr(jsonCASE)

cd('D:\Dropbox\OCD_DBS_JSON')
close all

switch jsonCASE
    case 1

        jsonfile = jsondecode(fileread('LeftHemisphere_afteresposures-1sttrial.json'));

    case 2

        jsonfile = jsondecode(fileread('LeftHemisphere_afterexposures-2ndtrial.json'));

    case 3

        jsonfile = jsondecode(fileread('LeftHemisphere_beforeexposures-1sttrial.json'));

    case 4

        jsonfile = jsondecode(fileread('LeftHemisphere_beforeexposures-2ndtrial.json'));

    case 5

        jsonfile = jsondecode(fileread('RightHemisphere_afterexposures-1sttrial.json'));

    case 6

        jsonfile = jsondecode(fileread('RightHemisphere_afterexposures-2ndtrial.json'));

    case 7

        jsonfile = jsondecode(fileread('RightHemisphere_beforeexposures-1sttrial.json'));

    case 8

        jsonfile = jsondecode(fileread('RightHemisphere_beforeexposures-2ndtrial.json'));

end


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
unfurl = reshape(smothLFPtrim,numel(smothLFPtrim),1);
normLFP = normalize(unfurl,'range');
repack = reshape(normLFP,6,72);
% plot(transpose(normSlfp))

% Reorder to max peak
averagePOWER = mean(repack,2);
[~ , high2low] = sort(averagePOWER , 'descend');

normSlfp = repack(high2low,:);

coloRS = [191, 15, 255;...
          194, 75, 210;...
          197, 135, 164;...
          200, 195, 119;...
          202, 225, 96;...
          203, 255, 73];
coloRSrgb = coloRS/255;

lineAlphas = linspace(0.9,0.2,6);
lineWidths = linspace(2,0.5,6);

normSlfpLt = transpose(normSlfp);

for lfpi = 1:width(normSlfpLt)

    hold on
    lalph = lineAlphas(lfpi);
    plot(tmpHzTRx,normSlfpLt(:,lfpi),'Color',[coloRSrgb(lfpi,:) lalph],...
        'LineWidth',lineWidths(lfpi))





end
hold off

senseELallS = senseELall(high2low,:);
senseELall2 = transpose(join(senseELallS,'-'));

xticks([0 35 70])
yticks([0 0.5 1])
ylabel('Normalized LFP magnitude')
xlabel('Frequency (Hz)')
legend(senseELall2)


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

















