%% Select subject JSON file

% NEED TO convert to A function
% TO DO:
% INPUTs:
% 1. JSON File location
% 2. Patient ID
% 3. Hemispheres collected [options: 'L','R','LR']
% 4. Number of runs [should hopefully always be 3]
% 5. double battery/single battery
% 6. Save data location

%%%%% CURRENT file runs the code on Patient 4

subjectFILE = 'pt5_0209_l_survey.json';

js = jsondecode(fileread(subjectFILE));

%%

% TO DO Create NEW CSV files with RUN / SIDE / CHANNEL 
% Old version only used the channel with MAX beta

% "1" = "L";
% "2" = "R";

channels = unique({js.LfpMontageTimeDomain.Channel}, 'stable'); % 

leng = numel({js.LfpMontageTimeDomain.Channel}); %
sides = {'LEFT', 'RIGHT'}; 
if any((contains(channels, sides{1}))) && any((contains(channels, sides{2}))) 
        doubleBattery = false; 
        %json has both left and right data 
        channelsLeft = channels(contains(channels, sides{1}));
        channelsRight = channels(contains(channels, sides{2}));
    else 
        doubleBattery = true; 
end 

% Time domain data
tdData = js.LfpMontageTimeDomain;
% All channels and runs LABELS
allcAr = {tdData.Channel};
allraw = {tdData.TimeDomainData};

% 4096 samples is the number of samples generated from computing the power
% spectral info from 250Hz signal - limited to 1-100 Hz.
numChans = 6;
numRuns = 3;
numSides = 1;
outDat = cell(numChans*numRuns*numSides,1);
outRUn = nan(numChans*numRuns*numSides,1);
outSide = cell(numChans*numRuns*numSides,1);
outChan = cell(numChans*numRuns*numSides,1);

combCount = 1;

% Loop throgh Lfp Montage Time Domain
if doubleBattery

    for ci = 1:length(channels)

        chanUnique = channels{ci};
        uniCruns = allraw(matches(allcAr,chanUnique));

        % Channel string info
        chanNameParts = split(chanUnique,'_');
        chanNAME = join(chanNameParts(1:3),'_');
        sideNAME = chanNameParts{4}(1);

        for ri = 1:length(uniCruns)

            uniRUN = uniCruns{ri};
            % Use known Fs = 250 Hz and limit to 0 - 100 Hz
            [pwR , freQ] = pspectrum(uniRUN, 250, 'FrequencyLimits', [0 100]);

            % betaBand = pwR(freQ > 12 & freQ < 33);
            outDat{combCount,1} = table(pwR,freQ,'VariableNames',{'Power','Frequency'});
            outRUn(combCount,1) = ri;
            outSide{combCount,1} = sideNAME;
            outChan{combCount,1} = chanNAME{1};

            combCount = combCount + 1;
            
        end
    end

else



end

outTABLE = table(outDat,outRUn,outSide,outChan,'VariableNames',...
    {'PF_Data','RunNum','SideID','ChanID'});

