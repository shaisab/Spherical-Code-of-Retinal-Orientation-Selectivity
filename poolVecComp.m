
function [pooledmap] = poolVecComp(pooledmap,vecCompAll,j)

    
%pooling and plotting vector fields in generalized (X,Y,Z) coordinates    
XYZvecCompAll=vecCompAll.XYZvecComp;
X=[];
Y=[];
Z=[];
VecX=[];
VecY=[];
VecZ=[];
angleXYZ=[];
for i=1:size(XYZvecCompAll,2)
X=[X;XYZvecCompAll(1,i).X'];
Y=[Y;XYZvecCompAll(1,i).Y'];
Z=[Z;XYZvecCompAll(1,i).Z'];
VecX=[VecX;XYZvecCompAll(1,i).VecX'];
VecY=[VecY;XYZvecCompAll(1,i).VecY'];
VecZ=[VecZ;XYZvecCompAll(1,i).VecZ'];
angleXYZ=[angleXYZ;XYZvecCompAll(1,i).angleXYZ'];
end




X=vecCompAll.XYZvecComp.X(j);
Y=vecCompAll.XYZvecComp.Y(j);
Z=vecCompAll.XYZvecComp.Z(j);
VecX=vecCompAll.XYZvecComp.VecX(j);
VecY=vecCompAll.XYZvecComp.VecY(j);
VecZ=vecCompAll.XYZvecComp.VecZ(j);
angleXYZ=vecCompAll.XYZvecComp.angleXYZ(j);



pooledmap.XYZvecComp.X=[];
pooledmap.XYZvecComp.Y=[];
pooledmap.XYZvecComp.Z=[];
pooledmap.XYZvecComp.VecX=[];
pooledmap.XYZvecComp.VecY=[];
pooledmap.XYZvecComp.VecZ=[];
pooledmap.XYZvecComp.angleXYZ=[];

pooledmap=struct;

pooledmap.XYZvecComp.X=[pooledmap.XYZvecComp.X;vecCompAll.XYZvecComp.X(j)];
pooledmap.XYZvecComp.Y=Y;
pooledmap.XYZvecComp.Z=Z;
pooledmap.XYZvecComp.VecX=VecX;
pooledmap.XYZvecComp.VecY=VecY;
pooledmap.XYZvecComp.VecZ=VecZ;
pooledmap.XYZvecComp.angleXYZ=angleXYZ;


%pooling and plotting vector fields in generalized STH coordinates
STHgenvecCompAll=vecCompAll.STHgenvecComp;
Xprime=[];
Yprime=[];
VecXprime=[];
VecYprime=[];
angleSTHGen=[];
for i=1:size(STHgenvecCompAll,2)
Xprime=[Xprime;STHgenvecCompAll(1,i).Xprime'];
Yprime=[Yprime;STHgenvecCompAll(1,i).Yprime'];
VecXprime=[VecXprime;STHgenvecCompAll(1,i).VecXprime'];
VecYprime=[VecYprime;STHgenvecCompAll(1,i).VecYprime'];
angleSTHGen=[angleSTHGen;STHgenvecCompAll(1,i).angleSTHGen'];
end

pooledmap.STHgenvecComp.Xprime=Xprime;
pooledmap.STHgenvecComp.Yprime=Yprime;
pooledmap.STHgenvecComp.VecXprime=VecXprime;
pooledmap.STHgenvecComp.VecYprime=VecYprime;
pooledmap.STHgenvecComp.angleSTHGen=angleSTHGen;


%pooling and plotting vector fields in STH coordinates
STHvecCompAll=vecCompAll.STHvecComp;
S=[];
TH=[];
VecS=[];
VecTH=[];
angleSTH=[];
for i=1:size(STHvecCompAll,2)
S=[S;STHvecCompAll(1,i).S'];
TH=[TH;STHvecCompAll(1,i).TH'];
VecS=[VecS;STHvecCompAll(1,i).VecS'];
VecTH=[VecTH;STHvecCompAll(1,i).VecTH'];
angleSTH=[angleSTH;STHvecCompAll(1,i).angleSTH'];
end

pooledmap.STHvecComp.S=S;
pooledmap.STHvecComp.TH=TH;
pooledmap.STHvecComp.VecS=VecS;
pooledmap.STHvecComp.VecTH=VecTH;
pooledmap.STHvecComp.angleSTH=angleSTH;


%pooling and plotting vector fields in UV coordinates
UVvecCompAll=vecCompAll.UVvecComp;
U=[];
V=[];
VecU=[];
VecV=[];
angleUV=[];
for i=1:size(UVvecCompAll,2)
U=[U;UVvecCompAll(1,i).U'];
V=[V;UVvecCompAll(1,i).V'];
VecU=[VecU;UVvecCompAll(1,i).VecU'];
VecV=[VecV;UVvecCompAll(1,i).VecV'];
angleUV=[angleUV;UVvecCompAll(1,i).angleUV'];
end

pooledmap.UVvecComp.U=U;
pooledmap.UVvecComp.V=V;
pooledmap.UVvecComp.VecU=VecU;
pooledmap.UVvecComp.VecV=VecV;
pooledmap.UVvecComp.angleUV=angleUV;


%pooling and plotting vector fields in UVstandard coordinates
UVStvecCompAll=vecCompAll.UVStvecComp;
USt=[];
VSt=[];
VecUSt=[];
VecVSt=[];
angleUVSt=[];
for i=1:size(UVStvecCompAll,2)
USt=[USt;UVStvecCompAll(1,i).USt'];
VSt=[VSt;UVStvecCompAll(1,i).VSt'];
VecUSt=[VecUSt;UVStvecCompAll(1,i).VecUSt'];
VecVSt=[VecVSt;UVStvecCompAll(1,i).VecVSt'];
angleUVSt=[angleUVSt;UVStvecCompAll(1,i).angleUVSt'];
end

pooledmap.UVStvecComp.USt=USt;
pooledmap.UVStvecComp.VSt=VSt;
pooledmap.UVStvecComp.VecUSt=VecUSt;
pooledmap.UVStvecComp.VecVSt=VecVSt;
pooledmap.UVStvecComp.angleUVSt=angleUVSt;

end


