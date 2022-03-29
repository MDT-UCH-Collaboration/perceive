function [] = ritalinLFPprocess_v2(analysis)

cd('C:\Users\John\Documents\GitHub\perceive\RitalinLFP_casereport')
preFILE = 'Report_Json_Session_Report_20220317T110442.json';
prejs = jsondecode(fileread(preFILE));


postFILE = 'Report_Json_Session_Report_20220317T114335.json';
postjs = jsondecode(fileread(postFILE));

close all
switch analysis

    case 'Bsense'

        darkCol1 = [1.0000 0.8235 0.1098];
        lightCol1 = [1.0000    0.9373    0.6902];
        darkCol2 = [0.4667    0.6745    0.1882];
        lightCol2 = [0.7529    0.9412    0.5059];

        [contactTpre] = processChanTab(prejs);
        preRing = matches(contactTpre.Type,'RING');
        [contactTpost] = processChanTab(postjs);
        postRing = matches(contactTpost.Type,'RING');

        tdDataPre = prejs.LfpMontageTimeDomain;
        tdDataPost = postjs.LfpMontageTimeDomain;

        % All channels and runs LABELS
        allcAr_Pre = {tdDataPre.Channel};
        allcAr_Pre1 = allcAr_Pre(preRing);
        allraw_Pre = {tdDataPre.TimeDomainData};
        allraw_Pre1 = allraw_Pre(preRing);

        allcAr_Post = {tdDataPost.Channel};
        allcAr_Post1 = allcAr_Post(postRing);
        allraw_Post = {tdDataPost.TimeDomainData};
        allraw_Post1 = allraw_Post(postRing);

        % 2a 2b 2c - Left ; 10a 10b 10c - Right
        % Plot Pre Left 2-1 and Post Left 2-1
        % Plot Pre Right 2-1 and Post Right 2-1

        allraw_PreL = allraw_Pre1(contains(allcAr_Pre1,'LEFT'));
        allcAr_PreL = allcAr_Pre1(contains(allcAr_Pre1,'LEFT'));

        maxPeak_PreL = findMaxBp(allraw_PreL);
        chanPreL = allcAr_PreL(maxPeak_PreL);
        rawPreL = allraw_PreL(maxPeak_PreL);

        allraw_PreR = allraw_Pre1(contains(allcAr_Pre1,'RIGHT'));
        allcAr_PreR = allcAr_Pre1(contains(allcAr_Pre1,'RIGHT'));

        % maxPeak_PreR = findMaxBp(allraw_PreR);
        chanPreR = allcAr_PreR(maxPeak_PreL);
        rawPreR = allraw_PreR(maxPeak_PreL);

        allraw_PostL = allraw_Post1(contains(allcAr_Post1,'LEFT'));
        allcAr_PostL = allcAr_Post1(contains(allcAr_Post1,'LEFT'));

        % maxPeak_PostL = findMaxBp(allraw_PostL);
        chanPostL = allcAr_PostL(maxPeak_PreL);
        rawPostL = allraw_PostL(maxPeak_PreL);

        allraw_PostR = allraw_Post1(contains(allcAr_Post1,'RIGHT'));
        allcAr_PostR = allcAr_Post1(contains(allcAr_Post1,'RIGHT'));

        chanPostR = allcAr_PostR(maxPeak_PreL);
        rawPostR = allraw_PostR(maxPeak_PreL);


        alldata = [rawPreL , rawPreR , rawPostL ,rawPostR];
        allChan = repmat(chanPreL,size(alldata));

        [allpwrsN3Pre , freQPre] = getNormPr(allChan , alldata);
        % [allpwrsN3Post , freQPost] = getNormPr(allcAr_Post1 , allraw_Post1);


        % plot(allpwrsN3Pre,'LineWidth',2.5)

        for pi = 1:4

            switch pi
                case 1 % Pre Left
                    plot(freQPre , allpwrsN3Pre(:,1) ,  'Color' , lightCol1, 'LineWidth',2.5)
                case 2 % Pre Right
                    plot(freQPre , allpwrsN3Pre(:,2) ,  'Color' , lightCol2, 'LineWidth',2.5)
                case 3 % Post Left
                    plot(freQPre , allpwrsN3Pre(:,3) ,  'Color' , darkCol1, 'LineWidth',2.5)
                case 4 % Post Right
                    plot(freQPre , allpwrsN3Pre(:,4) ,  'Color' , darkCol2, 'LineWidth',2.5)
            end
            hold on

        end



        xline(13,'Color','k','LineStyle','--')
        xline(30,'Color','k','LineStyle','--')
        title('ONE AND THREE')
        legend('+Ritalin LEFT','+Ritalin RIGHT', '-Ritalin LEFT', '-Ritalin RIGHT');
        xlabel('Frequency Hz')
        ylabel('Scaled power')
        axis square


    case 'stream'


        bsTD = {postjs.BrainSenseLfp.LfpData};

        allStreams = cell(length(bsTD),1);
        for bi = 1:length(bsTD)

            tmpStreamR = {bsTD{bi}.Right};
            tmpStreamL = {bsTD{bi}.Left};

            lLFP = zeros(length(tmpStreamL),1);
            rLFP = zeros(length(tmpStreamL),1);
            lSTIM = zeros(length(tmpStreamL),1);
            rSTIM = zeros(length(tmpStreamL),1);

            for rl = 1:2
                switch rl
                    case 1 % right
                        for bin = 1:length(tmpStreamR)
                            lLFP(bin) = tmpStreamR{bin}.LFP;
                            lSTIM(bin) = tmpStreamR{bin}.mA;
                        end
                    case 2 % left
                        for bin = 1:length(tmpStreamL)
                            rLFP(bin) = tmpStreamL{bin}.LFP;
                            rSTIM(bin) = tmpStreamL{bin}.mA;
                        end
                end
            end
            allStreams{bi} = table(lLFP,lSTIM,rLFP,rSTIM,'VariableNames',...
                {'LeftLFP','LeftmA','RightLFP','RightmA'});
        end


        % Clean up for publication
        rightDATA = allStreams{1};
        rightDATAa = rightDATA(41:height(rightDATA),:);

        figRIGHT = figure;
        left_color =  [0.4667    0.6745    0.1882];
        right_color = [0 0 0];
        set(figRIGHT,'defaultAxesColorOrder',[left_color; right_color]);

        xTimeSecs = 0:2:(height(rightDATAa)*2)-1; % seconds
        rightLfp = rightDATAa.RightLFP;
        yyaxis left

%         lfpsm = smoothdata(rightLfp,'movmean',15);

        plot(xTimeSecs,rightLfp,'LineWidth',2.5);
        ylabel('LFP Power')

        rightmA = rightDATAa.RightmA;
        yyaxis right
        plot(xTimeSecs,rightmA,'LineStyle','-.','LineWidth',2.5);
        ylabel('mA')

        xlabel('Time in seconds')
        title('Right STN')
        axis square

%         corrPLOT = figure;
%         plot(rightDATAa.RightLFP,rightDATAa.RightmA,'o')
%         hold on
        % Get coefficients of a line fit through the data.
%         coefficients = polyfit(rightDATAa.RightLFP,rightDATAa.RightmA, 1);
        % Create a new x axis with exactly 1000 points (or whatever you want).
%         xFit = linspace(min(rightDATAa.RightLFP), max(rightDATAa.RightLFP), 1000);
        % Get the estimated yFit value for each of those 1000 new x locations.
%         yFit = polyval(coefficients , xFit);
%         plot(xFit, yFit, 'r-', 'LineWidth', 2);
%         ylim([0 2.5])
        [rho , pval] = corr(rightDATAa.RightLFP,rightDATAa.RightmA);
%         title('r = -0.77 , p = 0.00001');
%         ylabel('mA')
%         xlabel('LFP power')


        


end


end




function [contactTable] = processChanTab(jsFILE)

channels = unique({jsFILE.LfpMontageTimeDomain.Channel}, 'stable');

strgroups = cellfun(@(x) strsplit(x, '_'), channels, 'UniformOutput', false);

contact1 = cell(length(strgroups),1);
contact2 = cell(length(strgroups),1);
contactT = cell(length(strgroups),1);
hemi = cell(length(strgroups),1);
for si = 1:length(strgroups)

    if length(strgroups{si}) == 5
        contact1{si} = strgroups{si}{1};
        contact2{si} = strgroups{si}{3};
        contactT{si} = strgroups{si}{5};
        hemi{si} = strgroups{si}{4};
    else
        contact1{si} = [strgroups{si}{1} strgroups{si}{2}];
        contact2{si} = [strgroups{si}{4} strgroups{si}{5}];
        contactT{si} = strgroups{si}{7};
        hemi{si} = strgroups{si}{6};
    end

end

contactTable = table(contact1,contact2, contactT, hemi, 'VariableNames',...
    {'Contact1','Contact2','Type','Hemi'});



end


function [allpwrsN3 , freQ] = getNormPr(channels , allraw)

allpwrs = zeros(4096,length(channels));
for ci = 1:length(channels)
    uniRUN = allraw{ci};
    [pwR , freQ] = pspectrum(uniRUN, 250, 'FrequencyLimits', [0.5 60]);
    pwRsm = smoothdata(pwR,"gaussian",300);
    allpwrs(:,ci) = pwRsm;
end

allpwrsN1 = allpwrs(:);
allpwrsN2 = normalize(allpwrsN1,'range');
allpwrsN3 = reshape(allpwrsN2, size(allpwrs));

end



function [maxPeakID] = findMaxBp(rawChannels)

allBetas = zeros(length(rawChannels),1);
maxPeakID = false(length(rawChannels),1);
for ci = 1:length(rawChannels)
    uniRUN = rawChannels{ci};
    [pwR , freQ] = pspectrum(uniRUN, 250, 'FrequencyLimits', [0.5 60]);
    betaB = mean(pwR(freQ > 13 & freQ < 30));
    allBetas(ci) = betaB;
end

[~ , maxID] = max(allBetas);
maxPeakID(maxID) = true;

end









