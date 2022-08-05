%%

clear
cd('C:\Users\John\Documents\GitHub\perceive')
%jsonFiles = 'Report_Json_Session_Report_20210604T094200.json';
jsonFiles = 'pt2_0604_rl_survey.json';
cd('C:\Users\John\Documents\GitHub\perceive')
js = jsondecode(fileread(jsonFiles));

channels = unique({js.LfpMontageTimeDomain.Channel}, 'stable');
leng = numel({js.LfpMontageTimeDomain.Channel});
sides = {'LEFT', 'RIGHT'};
startindex = 1;
if any((contains(channels, sides{1}))) && any((contains(channels, sides{2})))
    doubleBattery = false;
    %json has both left and right data
    channelsLeft = channels(contains(channels, sides{1}));
    channelsRight = channels(contains(channels, sides{2}));
    highestBetas = cell(2);
else
    doubleBattery = true;
end

close all
for  i = 1:2
    s1 = subplot(2, 1, i);
    % ax = gca;
    s1.FontSize = 14;

    %nexttile
    channels2 = channels(contains(channels,sides{i}));
    colors = 'krbgmc';
    maxValues = zeros(length(channels2), 1);
    for c=1:length(channels2)
        maxValues(c) = tempAvgPlot(startindex, leng, js, channels2{c}, colors(c));
        hold on;
    end

    if i == 1
        maxValuesLeft = maxValues;
    else
        maxValuesRight = maxValues;
    end

    [M, I] = max(maxValues);
    highestBeta = channels2{I};
    highestBetas{i} = highestBeta;

    %add all the plot info
    xline(13);
    xline(30);

    %formats the channels to be suitable for the legend
    legendChan = cell(6);
    oldchar = {'LEFT', 'RIGHT', '0', '1', '2', '3', '4', '5', '_', 'AND'};
    newchar = {''};
    oldnum = {'ZERO', 'ONE', 'TWO', 'THREE'};
    newnum = {'0', '1', '2', '3'};
    for b = 1:length(channels2)
        legendChan{b} = replace(channels2{b}, oldchar, newchar);
        legendChan{b} = replace(legendChan{b}, oldnum, newnum);
    end
    legend(legendChan{1}, legendChan{2}, legendChan{3}, legendChan{4}, legendChan{5}, legendChan{6}, 'FontSize', 10)


    title(["Frequency vs. Power", sides{i}])
    xlim([0 60])
    xlabel("Frequency (Hz)")
    ylabel("Power (mV)")

    disp(['The highest beta comes from contact pair ', legendChan{I}])
end

%%

close all
dir = readtable("PtDirectorySimple.xlsx");
mastertable = [];
%table with each run, contact, hemi and pt
for i = 1:height(dir)
    jsonFiles = dir.JsonSurvey{i};
    side = dir.Hemisphere{i};
    startIndex = dir.StartIndex(i);
    num = num2str(dir.Patient(i));
    hemi = upper(dir.Hemisphere(i));
    js = jsondecode(fileread(jsonFiles));
    outTABLE = extractSurveyData(js, side, startIndex);
    pt = repmat({append("pt", num)}, 18, 1);
    outTABLE.patient = pt;
    mastertable = [mastertable; outTABLE];
end

localtab = [];
y = zeros(4096,3);
mastertable2 = mastertable;
col = [];
%get only power data from mastertable
for x = 1:height(mastertable2)
    temp = mastertable.PF_Data{x};
    col = [col; temp.Power];
end
%normalize and repack
normcol = normalize(col, 'range'); %pack into 0 to 1
colm = reshape(col, 4096, height(mastertable2));

%put the normalized data into mastertable2
for x = 1:height(mastertable)
    mastertable2.PF_Data{x}.Power = colm(:,x);
end

for k = 1:height(mastertable2)
    %if we have the first run then we have that contact
    if mastertable2.RunNum(k) == 1
        t = mastertable2(k:k+2, :); %get all three runs for that contact
        for j = 1:3
            powerdata = t.PF_Data{j};
            y(:,j) = powerdata.Power;
            runav = mean(powerdata.Power);%single av for all power data
            runavmat(j) = runav; %store in matrix
        end
        ny1 = y(:); %takes y and makes it a column vector
        %ny2 = normalize(ny1, 'range'); %makes it 0 to 1
        ny3 = reshape(ny1, size(y)); %takes normalized data and reshapes to y
        ny3 = mean(ny3,2);
        runavall = mean(ny3); %average the average of all three runs, single number
        runst = std(ny3); %std dev of all three runs
        new = mastertable(k, :);
        new.PF_Data = [];
        new.RunNum = [];
        new.RunAvg = runavall;
        new.RunStd = runst;
        %new.avpwrdata = array2table(y);
        localtab = [localtab; new];

    end
end

%fix the ID categories for yticks
nameparts = cell(height(localtab), 3);
name = cell(height(localtab), 1);
for i = 1:height(localtab)
    oldchar = {'LEFT', 'RIGHT', '0', '1', '2', '3', '4', '5', '_', 'AND'};
    newchar = {''};
    oldnum = {'ZERO', 'ONE', 'TWO', 'THREE'};
    newnum = {'0', '1', '2', '3'};
    nameparts{i,1} = replace(localtab.ChanID{i}, oldchar, newchar);
    nameparts{i,1} = replace(nameparts{i,1}, oldnum, newnum);
    nameparts{i,1} = append((nameparts{i,1}(1)), '-', (nameparts{i,1}(2)));


    oldchar = {'LEFT', 'RIGHT'};
    newchar = {'L', 'R'};
    nameparts{i,2} = replace(localtab.SideID{i}, oldchar, newchar);

    nameparts{i,3} = localtab.patient{i};

    name{i,1} = append(nameparts{i,3}, ' ', nameparts{i,2}, ' ', nameparts{i,1});
end

ptnames = unique(cellstr(nameparts(:, 3)), 'stable');


avarray = localtab.RunAvg;
allmean = mean(avarray);
starray = localtab.RunStd;
allstd = mean(starray);
minarray = avarray - starray;
maxarray = avarray + starray;
%find which means exceed 1 std
avarray_more = avarray(avarray>(allmean+allstd));
indicies = find(avarray>(allmean+allstd));

y = [1:height(localtab); 1:height(localtab)];
x = transpose([avarray, maxarray]);
line(x,y, 'color', 'k', 'LineWidth', 0.5)
hold on
plot(avarray, y, 'Marker',"o", 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'LineStyle',"none")
hold on
plot(avarray_more, indicies, 'Marker',"o", 'MarkerFaceColor', 'red', 'LineStyle',"none")
xlim([-0.05 0.5])
xline(allmean+allstd, '-', {'1 standard'; 'deviation'}) %label this
%xline(allmean+2*allstd, '-', {'2sigma'})
xline(allmean, '-', {'Mean'})
xlabel("Normalized Power")
title("Average power range across all patients, hemispheres, and contacts")
text(0.4, 55, {'Blue Marker', 'Represents Average'})
sum(allmean+allstd<maxarray);
set(gca,'ytick',[1:height(localtab)],'yticklabel',name);
axis square

%%

sigruns = localtab(indicies,:)


y = zeros(4096,3);
t = tiledlayout(height(sigruns)+1, 1);

freq = mastertable.PF_Data{1}.Frequency;
for i = 1:height(sigruns)
    %copy on every loop to run rhough it
    mastercopy = mastertable;
    y = zeros(4096,3);
    %get the relevant PF data from sigurns
    mastercopy = mastercopy(contains(mastercopy.SideID, sigruns.SideID{i}) & contains(cellstr(mastercopy.patient), sigruns.patient{i}) & contains(mastercopy.ChanID, sigruns.ChanID{i}),:);
    %put into the y vector for storage
    y(:,1) = mastercopy.PF_Data{1}.Power;
    y(:,2) = mastercopy.PF_Data{2}.Power;
    y(:,3) = mastercopy.PF_Data{3}.Power;
    %get the average of the power for the relevant run
    y = mean(y,2);
    nY = normalize(y,'range');
    sNy = smoothdata(nY,'gaussian',30);
    %%%%PLOTTING
    nexttile
    plot(freq, sNy,'k','LineWidth',2)
    set(gca,'xtick',[])
    set(gca,'xtick',[])
    xlim([0 50])
    yticks([0, 0.5, 1])
    title(name(indicies(i)))

end

%PLOT NORMAL
y = zeros(4096,3);
%chose the first three runs of pt 1 left zero-three arbitrarily (it is
%indeed normal)
y(:,1) = mastertable.PF_Data{1}.Power;
y(:,2) = mastertable.PF_Data{2}.Power;
y(:,3) = mastertable.PF_Data{3}.Power;
y = mean(y,2);
nY = normalize(y,'range');
sNy = smoothdata(nY,'gaussian',30);
nexttile
plot(freq, sNy,'k','LineWidth',2)
xlim([0 50])
yticks([0, 0.5, 1])
title("Average Below 1 STD of Population Mean")