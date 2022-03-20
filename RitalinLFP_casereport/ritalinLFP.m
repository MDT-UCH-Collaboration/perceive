%% Read JSON File
cd('C:\Users\johna\Desktop\ritalin_LFP\PRE')
subjectFILE = 'Report_Json_Session_Report_20220317T110442.json';
js = jsondecode(fileread(subjectFILE));
%%
channels = unique({js.LfpMontageTimeDomain.Channel}, 'stable');

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

%%

% Time domain data
tdData = js.LfpMontageTimeDomain;

%%
% All channels and runs LABELS
allcAr = {tdData.Channel};
allraw = {tdData.TimeDomainData};

%%
close all
allpwrs = zeros(4096,length(channels));
for ci = 1:length(channels)
    uniRUN = allraw{ci};
    [pwR , freQ] = pspectrum(uniRUN, 250, 'FrequencyLimits', [0.5 60]);
    % Smooth
    pwRsm = smoothdata(pwR,"gaussian",300);
%     pwRsmN = normalize(pwRsm,'range');
%     plot(freQ,pwRsm)
%     hold on
    allpwrs(:,ci) = pwRsm;

end

%% Normalize
allpwrsN1 = allpwrs(:);
allpwrsN2 = normalize(allpwrsN1,'range');
%%
allpwrsN3 = reshape(allpwrsN2, size(allpwrs));
plot(allpwrsN3)
