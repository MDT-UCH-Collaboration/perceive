% monoTab = js.Impedance.Hemisphere.SessionImpedance.Monopolar;
% 
% % Electrode2 = get number that was used for stimulation
% % ohn, could we double check the JSON file for this to make sure we're right? It's subject #6
% % PW, Frequency, amplitude
% 
% finGroup = js.Groups.Final;
% settINGS = finGroup.ProgramSettings;
% ratEE = settINGS.RateInHertz;
% sensCCg = settINGS.SensingChannel;
% pulseWt = sensCCg.PulseWidthInMicroSecond;
% 
% % ResultValue
% imped = monoTab(contains({monoTab.Electrode2},'ElectrodeDef.FourElectrodes_0'),:);
% impedU = imped.ResultValue;
% 
% averageAmp = 2.1;
% 
% teed = averageAmp^2 * impedU * ratEE * pulseWt;
% % 76776336 = 2.1^2 * 1488 * 90 * 130


%%

computerName = getenv("COMPUTERNAME");

switch computerName
    case 'DESKTOP-EGFQKAI'

        cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020')
        load('PatParms.mat','patParams')
        mainDIR = 'C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';

    otherwise

        cd('E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020')
        load('PatParms.mat','patParams')
        mainDIR = 'E:\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020';


end
userDIRs = [mainDIR,'\Data\SPPD'];
userDIRe = '\JSON_LFP';



stimmA = [3 , 3.2 , 2.0 , 2.6 , 3.1 , 2.0 , 4.4, 3.8, 3.4, 4.6, 6, 2.6, 2.0, 2.5, 1.0, 1.8, 1.8, 1.5 ];
elecID = [1 1 1 2 2 1 2 1 1 2 2 2 1 1 1 2 2 2];

pifields = fieldnames(patParams);
piCOUNT = 1;
allTEED = zeros(length(stimmA),1);

for pi = 1:length(pifields)

    tmpPI = patParams.(pifields{pi});
    tmpPIname = pifields{pi}(2:end);

    hemifields = fieldnames(tmpPI);

    for hi = 1:length(hemifields)
        if matches(hemifields{hi},'R')
            hemiLOC = '\Right';
        else
            hemiLOC = '\Left';
        end
        tmpPIloc = [userDIRs,tmpPIname,userDIRe,hemiLOC];
        cd(tmpPIloc)

        jsName = patParams.(pifields{pi}).(hemifields{hi}).json;
        js = jsondecode(fileread(jsName));

        finGroup = js.Groups.Final;

        activeGROUP = find([finGroup.ActiveGroup]);
        settINGS = finGroup(activeGROUP).ProgramSettings;
        if isfield(settINGS,'RateInHertz')
            ratEE = settINGS.RateInHertz;
        else
            ratEE = settINGS.SensingChannel.RateInHertz;
        end
        sensCCg = settINGS.SensingChannel;
        pulseWt = sensCCg.PulseWidthInMicroSecond;

        testIMPD = js.Impedance.Hemisphere;
        if isempty(fieldnames(testIMPD))
            allTEED(piCOUNT) = nan;
            piCOUNT = piCOUNT + 1;
        else

            if length([js.Impedance.Hemisphere]) ~= 1

                if matches(hemifields{hi},'R')
                    rowID = contains({js.Impedance.Hemisphere.Hemisphere},'Right');
                    monoTab = js.Impedance.Hemisphere(rowID).SessionImpedance.Monopolar;
                else
                    rowID = contains({js.Impedance.Hemisphere.Hemisphere},'Left');
                    monoTab = js.Impedance.Hemisphere(rowID).SessionImpedance.Monopolar;
                end
            else
                monoTab = js.Impedance.Hemisphere.SessionImpedance.Monopolar;
            end

            eleCtest = contains(extractAfter({monoTab.Electrode2},'_'),num2str(elecID(piCOUNT)));

            imped = monoTab(eleCtest,:);
            impedU = imped.ResultValue;

            averageAmp = stimmA(piCOUNT);

            allTEED(piCOUNT) = averageAmp^2 * impedU * ratEE * pulseWt;

            piCOUNT = piCOUNT + 1;
        end

    end

end



