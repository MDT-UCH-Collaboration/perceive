% Pring al pdfs


caseS = {'1','2','3','4','5','6','7','7','8','9','9','10'};
hemiS = {'R','L','L','R','R','R','L','R','L','R','L','R'};
evFLAG = [1 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 1 , 1 , 1];
savePDFloc = 'C:\Users\John\Desktop\pdfSsave';

for ci = 1:length(caseS)

    makeAlexPerceptFigures_v3(1,caseS{ci},hemiS{ci},evFLAG(ci))
    cd(savePDFloc)
    exportgraphics(gcf,[caseS{ci},'_',hemiS{ci},'.pdf'])

end