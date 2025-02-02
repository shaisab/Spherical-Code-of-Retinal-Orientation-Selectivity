function []= VecComparisonAll(cellType, fieldType)

% cellType can be either 'ONDS', 'ONOFFDS', 'OFFDS', or 'retroONDS'. 

%VecComparison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script compares the vector fields obtained from the elastic model
%   with that from Shai's data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input vector field data from model.
S1=dlmread(['Sec1VecField_',fieldType]);
S2=dlmread(['Sec2VecField_',fieldType]);
S3=dlmread(['Sec3VecField_',fieldType]);
S4=dlmread(['Sec4VecField_',fieldType]);


%number of points discretized over
n=size(S1,1);

%Sector 1 Vector Field
U1=S1(1:n,1:n);
V1=S1(1:n,n+1:2*n);
VecU1=S1(1:n,2*n+1:3*n);
VecV1=S1(1:n,3*n+1:4*n);

%Sector 2 Vector Field
U2=S2(1:n,1:n);
V2=S2(1:n,n+1:2*n);
VecU2=S2(1:n,2*n+1:3*n);
VecV2=S2(1:n,3*n+1:4*n);

%Sector 3 Vector Field
U3=S3(1:n,1:n);
V3=S3(1:n,n+1:2*n);
VecU3=S3(1:n,2*n+1:3*n);
VecV3=S3(1:n,3*n+1:4*n);

%Leaf 4 Vector Field
U4=S4(1:n,1:n);
V4=S4(1:n,n+1:2*n);
VecU4=S4(1:n,2*n+1:3*n);
VecV4=S4(1:n,3*n+1:4*n);

% DirectionSelective Cells
[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

namesMap=getFileList(filedir,'map_',0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
    load(num2str(namesMap{i}));
    end
end

%load('map_May_23_2014.mat');

data=map.SzMapOutputDS;

U=data(:,1);
V=data(:,2);
VecU=data(:,3);
VecV=data(:,4);
Sec=data(:,5);

switch cellType
    case 'ONDS' 
type=map.ONDS;
    case 'ONOFFDS' 
type=map.ONOFFDS;
    case 'OFFDS' 
type=map.OFFDS;
    case 'retroONDS' 
type=map.retroONDS;
end


n=size(U,1);

figure;
hold on;
j=0;
for i=1:n,

    if type(i)==1
        j=j+1;   
   if Sec(i)==1,
       
       temp=(U(i)-U1).^2+(V(i)-V1).^2;
       [ind1,ind2] =find(temp==min(min(temp)));
       a=ind1(1);
       b=ind2(1);
       
       
       angleDS(j)=acos((VecU(i)*VecU1(a,b)+VecV(i)*VecV1(a,b))/...
           (sqrt(VecU(i)^2+VecV(i)^2)*(sqrt(VecU1(a,b)^2+VecV1(a,b)^2))));
       
       Ucomp(i)=U1(a,b);
       Vcomp(i)=V1(a,b);
       VecUcomp(i)=VecU1(a,b);
       VecVcomp(i)=VecV1(a,b);
       
   elseif Sec(i)==2,
       
       temp=(U(i)-U2).^2+(V(i)-V2).^2;
       [ind1,ind2] =find(temp==min(min(temp)));
       a=ind1(1);
       b=ind2(1);
       
       
       angleDS(j)=acos((VecU(i)*VecU2(a,b)+VecV(i)*VecV2(a,b))/...
           (sqrt(VecU(i)^2+VecV(i)^2)*(sqrt(VecU2(a,b)^2+VecV2(a,b)^2))));
       
       Ucomp(i)=U2(a,b);
       Vcomp(i)=V2(a,b);
       VecUcomp(i)=VecU2(a,b);
       VecVcomp(i)=VecV2(a,b);
       
       
   elseif Sec(i)==3,
       temp=(U(i)-U3).^2+(V(i)-V3).^2;
       [ind1,ind2] =find(temp==min(min(temp)));
       a=ind1(1);
       b=ind2(1);
       
       
       angleDS(j)=acos((VecU(i)*VecU3(a,b)+VecV(i)*VecV3(a,b))/...
           (sqrt(VecU(i)^2+VecV(i)^2)*(sqrt(VecU3(a,b)^2+VecV3(a,b)^2))));
       
       Ucomp(i)=U3(a,b);
       Vcomp(i)=V3(a,b);
       VecUcomp(i)=VecU3(a,b);
       VecVcomp(i)=VecV3(a,b);
   else
       temp=(U(i)-U4).^2+(V(i)-V4).^2;
       [ind1,ind2] =find(temp==min(min(temp)));
       a=ind1(1);
       b=ind2(1);
       
       
       angleDS(j)=acos((VecU(i)*VecU4(a,b)+VecV(i)*VecV4(a,b))/...
           (sqrt(VecU(i)^2+VecV(i)^2)*(sqrt(VecU4(a,b)^2+VecV4(a,b)^2))));
       
       Ucomp(i)=U4(a,b);
       Vcomp(i)=V4(a,b);
       VecUcomp(i)=VecU4(a,b);
       VecVcomp(i)=VecV4(a,b);
       
   end
end
end
dlmwrite(['angleDifDS_',fieldType],angleDS);
scale=0.2;
quiver(Ucomp,Vcomp,VecUcomp,VecVcomp);
quiver(U(type),V(type),VecU(type),VecV(type),scale);
hold off;
end

