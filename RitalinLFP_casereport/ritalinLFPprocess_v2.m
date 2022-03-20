function [] = ritalinLFPprocess_v2()

cd('C:\Users\johna\Desktop\ritalin_LFP\PRE')
preFILE = 'Report_Json_Session_Report_20220317T110442.json';
prejs = jsondecode(fileread(preFILE));

cd('C:\Users\johna\Desktop\ritalin_LFP\POST')
postFILE = 'Report_Json_Session_Report_20220317T114335.json';
postjs = jsondecode(fileread(postFILE));

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
allcAr_PreL2 = allraw_PreL(contains(allcAr_PreL,'TWO'));

allraw_PreR = allraw_Pre1(contains(allcAr_Pre1,'RIGHT'));
allcAr_PreR = allcAr_Pre1(contains(allcAr_Pre1,'RIGHT'));
allcAr_PreR2 = allraw_PreR(contains(allcAr_PreR,'TWO'));

allraw_PostL = allraw_Post1(contains(allcAr_Post1,'LEFT'));
allcAr_PostL = allcAr_Post1(contains(allcAr_Post1,'LEFT'));
allcAr_PostL2 = allraw_PostL(contains(allcAr_PostL,'TWO'));

allraw_PostR = allraw_Post1(contains(allcAr_Post1,'RIGHT'));
allcAr_PostR = allcAr_Post1(contains(allcAr_Post1,'RIGHT'));
allcAr_PostR2 = allraw_PostR(contains(allcAr_PostR,'RIGHT'));


[allpwrsN3Pre , freQPre] = getNormPr(allcAr_Pre1 , allraw_Pre1);
[allpwrsN3Post , freQPost] = getNormPr(allcAr_Post1 , allraw_Post1);


plot(allpwrsN3Pre)
hold    
plot(allpwrsN3Post)


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


