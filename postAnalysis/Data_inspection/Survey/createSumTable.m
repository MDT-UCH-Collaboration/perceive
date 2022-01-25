function [sumTABLEp, sumTABLEh] = createSumTable(outTABLE)
index = [1,4,7,10,13,16; 2,5,8,11,14,17; 3,6,9,12,15,18];
M = zeros(1,6); 
for j = 1:length(index)

    runs = cell(3,1);
    n = 1;
    means = zeros(860, 3); 

    for i = 1:3
        x = index(i,j);
        powerdata = outTABLE.PF_Data{x};
        fbeta = powerdata.Frequency > 12 & powerdata.Frequency < 33;
        %get raw beta data and put it in a table 
        betapwr = powerdata.Power(fbeta);
        means(:, i) = betapwr; 
        %modify data for KS test 
        betapwr = abs(log10(betapwr));
        betapwr = downsample(betapwr, 10);
        runs{n} = betapwr;
        n=n+1;
    end

    means = mean(means, 1); 
    M(j) = mean(means, 2);
    %chan{j} = outTABLE.ChanID(index(1,j)); 
    %P12 is the p value for the KS difference between runs 1 and 2
    %H12 is a logical value in which 0 is non sig diff and 1 is non sig diff
    [H12, P12] = kstest2(runs{1}, runs{2});
    [H13, P13] = kstest2(runs{1}, runs{3});
    [H23, P23] = kstest2(runs{2}, runs{3});

    %create summary table with P values
    sumTABLEptemp = array2table([M(j); P12; P13; P23], "VariableNames", outTABLE.ChanID(index(1,j))); %holds current values for channel
    if j == 1
        sumTABLEp = sumTABLEptemp; %in the first channel, the first channel is the whole table
    else
        sumTABLEp = [sumTABLEp sumTABLEptemp]; %for all other channels, add the new channels to the old ones
    end
    %add the mean beta variable 
    
% 
    %create summary table with H values
    sumTABLEhtemp = array2table([M(j); H12; H13; H23], "VariableNames", outTABLE.ChanID(index(1,j)));
    if j == 1
        sumTABLEh = sumTABLEhtemp;
    else
        sumTABLEh = [sumTABLEh sumTABLEhtemp];
    end 
end
sumTABLEp.Properties.RowNames = ["Max Beta", "Run 1 vs 2", "Run 1 vs 3", "Run 2 vs 3"];
sumTABLEh.Properties.RowNames = ["Max Beta", "Run 1 vs 2", "Run 1 vs 3", "Run 2 vs 3"]; 

[~, y] = sort(table2array(sumTABLEp(1, :)), "descend"); %x is values in correct order %y is indicies %use tilde to suppress specific output
sumTABLEp = sumTABLEp(:,y); 
sumTABLEh = sumTABLEh(:,y); 

end 