function [sumTABLEpmaster, sumTABLEhmaster] = createSumTable(outTABLE)
%used to go through all the data 
index = [1,4,7,10,13,16; 2,5,8,11,14,17; 3,6,9,12,15,18];
edges = [5 9; 9 12; 12 33; 60 90];  %frequency thresholds
dwnsample = [2; 2; 10; 14]; %gets raw data down to ~90 samples
bands = {'theta','alpha','beta','gamma'};

for k = 1:length(edges)
    M = zeros(1,6);
    
    for j = 1:length(index) %this loop runs through all the channels containing data
        
        runs = cell(3,1); %stores the avg power
        means = []; %stores all the raw data for beta to be avgeraged
        n = 1;
        
        %this loop runs through each channel data in outTABLE
        %impt to note in outTABLE the channels are grouped together
        %in the struct, not the case
        for i = 1:3 %moves across row in index to get data (1,2,3) , (4,5,6)
            x = index(i,j);
            
            %get beta (12-33hz)
            powerdata = outTABLE.PF_Data{x};
            %find freq vales of interest
            fbeta = powerdata.Frequency > edges(k,1) & powerdata.Frequency < edges(k,2);
            %ftheta = powerdata.Frequency > 5 & powerdata.Frequency < 9;
            %get raw beta data and put it in a table based on freq vals
            
            betapwr = powerdata.Power(fbeta);
            %thetapwr = powerdata.Power(ftheta);
            
            means(:, i) = betapwr;
            %meanst(:, i) = thetapwr;
            %modify data for KS test
            
            betapwr = abs(log10(betapwr));
            %thetapwr = abs(log10(thetapwr));
            
            %downsample for KS test not for means calculation
            betapwr = downsample(betapwr, 10);
            %thetapwr = downsample(thetapwr, 2);
            
            runs{n} = betapwr;
            %runst{n} = thetapwr;
            
            n=n+1;
        end
        
        means = mean(means, 1); %sums down column
        M(j) = mean(means, 2); %sums across 1 to get 1 number
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
        sumTABLEp.Properties.RowNames = [append("Max ", bands{k}), append("r1v2 ", bands{k}), append("r1v3 ", bands{k}), append("r2v3 ", bands{k})];
        sumTABLEh.Properties.RowNames = [append("Max ", bands{k}), append("r1v2 ", bands{k}), append("r1v3 ", bands{k}), append("r2v3 ", bands{k})];
    if k == 1 %save the table just created in master variable 
        sumTABLEpmaster = sumTABLEp; 
        sumTABLEhmaster = sumTABLEh; 
    else %append w new table if we can 
        %sumTABLEpmaster = [sumTABLEpmaster; sumTABLEp]; 
        sumTABLEpmaster = [sumTABLEpmaster; sumTABLEp]; 
        sumTABLEhmaster = [sumTABLEhmaster; sumTABLEh]; 
    end 
       
end

% [~, y] = sort(table2array(sumTABLEp(1, :)), "descend"); %x is values in correct order %y is indicies %use tilde to suppress specific output
% sumTABLEp = sumTABLEp(:,y);
% sumTABLEh = sumTABLEh(:,y);

end