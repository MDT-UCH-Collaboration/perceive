load('jattable_allbeta.mat','finalTable')
% patient
patLIST = string(finalTable.PatID);
patID = unique(patLIST);
figtab = table;
i = 0;
for pi = 1:length(patID)
    patTable = finalTable(matches(patLIST,patID(pi)),:);
    % side of brain
    sideList = string(patTable.SideID);
    sideID = unique(sideList);
    for si = 1:length(sideID)
        i = i + 1;
        sideTable = patTable(matches(sideList,sideID(si)),:);
        %find the highest two beta and store them in a new table
        betasort = sortrows(sideTable, 'BetaPw', 'descend');
        twoBeta = betasort(1:2, :);
        betaDiff = twoBeta.BetaHz(1) - twoBeta.BetaHz(2);
        figtab.PatID{i} = betasort.PatID{1};
        figtab.SideID{i} = betasort.SideID{1};
        figtab.BetaDiff(i) = betaDiff;
        figtab.bigbeta(i) = twoBeta.BetaHz(1);
        figtab.smallbeta(i) = twoBeta.BetaHz(2);
        contList = string(sideTable.Hemi);
        % Find highest and second highest beta power
        % Subject Peak Hz between the two contacts
    end
end
%creating the figure
%create x axis categories
for j = 1:height(figtab)
    ptside(j) = append(figtab.PatID{j}, ' ', extractBetween(figtab.SideID{j}, 1,1));
end

%%

figure
xdata= categorical(ptside)';
x = linspace(1, 9, 9);
ydata = (figtab.BetaDiff);
% plot(xdata, ydata, 'Marker', "o", 'LineStyle',"none", 'MarkerSize', 9, 'MarkerFaceColor', "red", 'MarkerEdgeColor', "red")
% yline(0, 'LineStyle', "--")
% ylabel("Frequency Difference")
% xlabel("Patient and Hemisphere")
figure
plot(figtab.bigbeta, x, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'Marker', "o", 'LineStyle',"none", 'MarkerSize', 9)
hold on
plot(figtab.smallbeta, x, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue', 'Marker', "o", 'LineStyle',"none", 'MarkerSize', 9)
xlabel("Frequency (Hz)")
ylabel("Patient and Hemisphere")
%legend("Hz of MBC", "Hz of 2nd MBC", "Mean Frequency of MBC", "Mean Frequency of 2nd MBC")

%%

for k = 1:height(figtab)
    p1 = [figtab.bigbeta(k), figtab.smallbeta(k)];
    p2 = [x(k), x(k)];
    hold on
    plot(p1, p2, 'Color', 'black')
end

%%
yticklabels(xdata)
ylim([0 10])
bigbetamean = mean(figtab.bigbeta);
smallbetamean = mean(figtab.smallbeta);
xline(bigbetamean, 'Color', 'red', 'LineStyle', '--');
xline(smallbetamean, 'color', 'blue', 'LineStyle', '--');
% legend("Hz of MBC", "Hz of 2nd MBC", "Mean Frequency of MBC", "Mean Frequency of 2nd MBC")


axis square
