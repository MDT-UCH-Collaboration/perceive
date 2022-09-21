
load('jattable_allbeta.mat','finalTable')

% patient
patLIST = string(finalTable.PatID);
patID = unique(patLIST);

for pi = 1:length(patID)

    patTable = finalTable(matches(patLIST,patID(pi)),:);

    % side of brain
    sideList = string(patTable.SideID);
    sideID = unique(sideList);

    for si = 1:length(sideID)
        sideTable = patTable(matches(sideList,sideID(si)),:);

        contList = string(sideTable.Hemi);

        % Find highest and second highest beta power

        % Subject Peak Hz between the two contacts



    end
end