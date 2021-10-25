
mainLOC = 'C:\Users\johna\Dropbox\Publications_Meta\InProgress\ABaumgartner_Percept2020\Data\';
subID = 2;
subNUM = ['SPPD' , num2str(subID)];
actLOC = [mainLOC,subNUM,'\ACT_data'];

cd(actLOC)

csvTABLE = csvFIND()




































function [outCELL] = csvFIND()

csvDir = dir('*.csv');
csvNames = {csvDir.name};

fid = fopen(csvNames{1});
tline = fgetl(fid);
lineCount = 1;
outCELL = cell(1,1000000);
while ischar(tline)
    outCELL{lineCount} = tline;
    lineCount = lineCount + 1;
    tline = fgetl(fid);
end
fclose(fid);

outCELL = outCELL(cellfun(@(x) ~isempty(x) , outCELL , 'UniformOutput',false));

end

