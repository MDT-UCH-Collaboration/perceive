[lvnaBound,lvnaVol] = boundary(lvna,0);

hold on
trisurf(lvnaBound,lvna(:,1),lvna(:,2),lvna(:,3),'Facecolor','red','FaceAlpha',0.1,'EdgeColor','none')

% drawMesh(lvna,k)



hold on

[lstnBound,lstnVol] = boundary(lstn,0);

hold on
trisurf(lstnBound,lstn(:,1),lstn(:,2),lstn(:,3),'Facecolor','blue','FaceAlpha',0.1,'EdgeColor','none')

% drawMesh(lstn,k)
title([num2str(lvnaVol) , ' ', num2str(lstnVol)])