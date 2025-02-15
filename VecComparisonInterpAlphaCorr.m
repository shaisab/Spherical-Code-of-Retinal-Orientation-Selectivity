function []= VecComparisonInterpAlphaCorr(discPoints,fieldType,OSIlow,OSIhigh)

% cellType can be either 'ONDS', 'ONOFFDS', 'OFFDS', or 'retroONDS'.
% fieldType can be either 'ASC', 'PSC', 'LSC', 'cardinalUp', etc.
%VecComparison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script compares the vector fields obtained from the elastic model
%   with that from Shai's data. We compare the data in the four releveant
%   coordinate systems for the problem. We write the data to four seperate
%   files called
%`
%   1. angleDifUV
%   2. angleDifSTH
%   3. angleDifGSTH
%   4. angleDifXYZ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Physical Parameters.
% This data is imported from data exported from the mapping script. These
% values correspond to measured values on the flattened retina.

phys=dlmread('physAlphaCorr');

R=phys(1);
M=phys(2); %This should be normalized to 1.

m1=phys(3);
m2=phys(4);
m3=phys(5);
m4=phys(6);

nu=phys(7);

a1=phys(8);
a2=phys(9);
a3=phys(10);
a4=phys(11);

%% Input vector field data from model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Here we input the data from the model. Speficially, we need the data
%   generated by the mapping file and the vector field scripts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    Sector1=dlmread(['Sec1VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType]);
    Sector2=dlmread(['Sec2VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType]);
    Sector3=dlmread(['Sec3VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType]);
    Sector4=dlmread(['Sec4VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType]);
catch
    disp('no such vector field exist');
    return
end

%number of points discretized over
n=size(Sector1,1);

%% Sector 1 Vector Field Data
U1=Sector1(1:n,1:n);
V1=Sector1(1:n,n+1:2*n);
VecU1=Sector1(1:n,2*n+1:3*n);
VecV1=Sector1(1:n,3*n+1:4*n);

RHO1=Sector1(1:n,4*n+1:5*n);
F1=Sector1(1:n,5*n+1:6*n);
RHOS1=Sector1(1:n,6*n+1:7*n);
RHOTH1=Sector1(1:n,7*n+1:8*n);
FS1=Sector1(1:n,8*n+1:9*n);
FTH1=Sector1(1:n,9*n+1:10*n);
S1=Sector1(1:n,10*n+1:11*n);
TH1=Sector1(1:n,11*n+1:12*n);
VecS1=Sector1(1:n,12*n+1:13*n);
VecTH1=Sector1(1:n,13*n+1:14*n);

Xprime1=Sector1(1:n,14*n+1:15*n);
Yprime1=Sector1(1:n,15*n+1:16*n);
VecXprime1=Sector1(1:n,16*n+1:17*n);
VecYprime1=Sector1(1:n,17*n+1:18*n);

X1=Sector1(1:n,18*n+1:19*n);
Y1=Sector1(1:n,19*n+1:20*n);
Z1=Sector1(1:n,20*n+1:21*n);
VecX1=Sector1(1:n,21*n+1:22*n);
VecY1=Sector1(1:n,22*n+1:23*n);
VecZ1=Sector1(1:n,23*n+1:24*n);

%% Sector 2 Vector Field Data
U2=Sector2(1:n,1:n);
V2=Sector2(1:n,n+1:2*n);
VecU2=Sector2(1:n,2*n+1:3*n);
VecV2=Sector2(1:n,3*n+1:4*n);

RHO2=Sector2(1:n,4*n+1:5*n);
F2=Sector2(1:n,5*n+1:6*n);
RHOS2=Sector2(1:n,6*n+1:7*n);
RHOTH2=Sector2(1:n,7*n+1:8*n);
FS2=Sector2(1:n,8*n+1:9*n);
FTH2=Sector2(1:n,9*n+1:10*n);
S2=Sector2(1:n,10*n+1:11*n);
TH2=Sector2(1:n,11*n+1:12*n);
VecS2=Sector2(1:n,12*n+1:13*n);
VecTH2=Sector2(1:n,13*n+1:14*n);

Xprime2=Sector2(1:n,14*n+1:15*n);
Yprime2=Sector2(1:n,15*n+1:16*n);
VecXprime2=Sector2(1:n,16*n+1:17*n);
VecYprime2=Sector2(1:n,17*n+1:18*n);

X2=Sector2(1:n,18*n+1:19*n);
Y2=Sector2(1:n,19*n+1:20*n);
Z2=Sector2(1:n,20*n+1:21*n);
VecX2=Sector2(1:n,21*n+1:22*n);
VecY2=Sector2(1:n,22*n+1:23*n);
VecZ2=Sector2(1:n,23*n+1:24*n);

%% Sector 3 Vector Field Data
U3=Sector3(1:n,1:n);
V3=Sector3(1:n,n+1:2*n);
VecU3=Sector3(1:n,2*n+1:3*n);
VecV3=Sector3(1:n,3*n+1:4*n);

RHO3=Sector3(1:n,4*n+1:5*n);
F3=Sector3(1:n,5*n+1:6*n);
RHOS3=Sector3(1:n,6*n+1:7*n);
RHOTH3=Sector3(1:n,7*n+1:8*n);
FS3=Sector3(1:n,8*n+1:9*n);
FTH3=Sector3(1:n,9*n+1:10*n);
S3=Sector3(1:n,10*n+1:11*n);
TH3=Sector3(1:n,11*n+1:12*n);
VecS3=Sector3(1:n,12*n+1:13*n);
VecTH3=Sector3(1:n,13*n+1:14*n);

Xprime3=Sector3(1:n,14*n+1:15*n);
Yprime3=Sector3(1:n,15*n+1:16*n);
VecXprime3=Sector3(1:n,16*n+1:17*n);
VecYprime3=Sector3(1:n,17*n+1:18*n);

X3=Sector3(1:n,18*n+1:19*n);
Y3=Sector3(1:n,19*n+1:20*n);
Z3=Sector3(1:n,20*n+1:21*n);
VecX3=Sector3(1:n,21*n+1:22*n);
VecY3=Sector3(1:n,22*n+1:23*n);
VecZ3=Sector3(1:n,23*n+1:24*n);

%% Sector 4 Vector Field Data
U4=Sector4(1:n,1:n);
V4=Sector4(1:n,n+1:2*n);
VecU4=Sector4(1:n,2*n+1:3*n);
VecV4=Sector4(1:n,3*n+1:4*n);

RHO4=Sector4(1:n,4*n+1:5*n);
F4=Sector4(1:n,5*n+1:6*n);
RHOS4=Sector4(1:n,6*n+1:7*n);
RHOTH4=Sector4(1:n,7*n+1:8*n);
FS4=Sector4(1:n,8*n+1:9*n);
FTH4=Sector4(1:n,9*n+1:10*n);
S4=Sector4(1:n,10*n+1:11*n);
TH4=Sector4(1:n,11*n+1:12*n);
VecS4=Sector4(1:n,12*n+1:13*n);
VecTH4=Sector4(1:n,13*n+1:14*n);

Xprime4=Sector4(1:n,14*n+1:15*n);
Yprime4=Sector4(1:n,15*n+1:16*n);
VecXprime4=Sector4(1:n,16*n+1:17*n);
VecYprime4=Sector4(1:n,17*n+1:18*n);

X4=Sector4(1:n,18*n+1:19*n);
Y4=Sector4(1:n,19*n+1:20*n);
Z4=Sector4(1:n,20*n+1:21*n);
VecX4=Sector4(1:n,21*n+1:22*n);
VecY4=Sector4(1:n,22*n+1:23*n);
VecZ4=Sector4(1:n,23*n+1:24*n);

%% Collection from Shai's data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   We now input measured Data from the retina in the flattened
%   coordinates.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [filedir,go]=bbdirselector('select data folder',cd);
% if ~go disp('no folder selected'); return; end
filedir=cd;

namesMap=getFileList(filedir,'mapAlphaCorr_',0,'anywhere');
for i=1:size(namesMap,2)
    if strfind(namesMap{i},'.mat')
        load(num2str(namesMap{i}));
    end
end

%load('map_May_27_2014_2.mat');

%  data=map.SzMapOutputDS; % for all imaged cells and DS
data=map.SzMapOutputOS_bars_OFF; % for OS

Utemp=data(:,1);
Vtemp=data(:,2);
VecUtemp=data(:,3);
VecVtemp=data(:,4);
Sec=data(:,5);


type=map.isOS .* (map.OSI>OSIlow) .* (map.OSI<OSIhigh);      % include all OS cells with OSI between OSI low and OSI high in the analysis
% type=ones(size(map.isOS,1),1);  % for all imaged cells

n=size(Utemp,1);

injectionSite={};
i=0;

%% Locating closest point on numerical mesh, converting to other vector fields and calculating angular difference.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   We now locate the point on the discrete mesh that is closest to the
%   point in Shai's data. We can take the closest point in the l^2 norm. We
%   record the index of the data points in the variables a(i), b(i) which
%   we can use to find the cooresponding point in other coordinate
%   systems. We then convert the data to other coordinate systems and
%   record the angle between the measured data and the hypothesized model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

% figure;
% hold on;
% surf(U1,V1,zeros(size(U1,1),size(U1,1)), 'EdgeColor', 'none');
% surf(U2,V2,zeros(size(U2,1),size(U2,1)), 'EdgeColor', 'none');
% surf(U3,V3,zeros(size(U3,1),size(U3,1)), 'EdgeColor', 'none');
% surf(U4,V4,zeros(size(U4,1),size(U4,1)), 'EdgeColor', 'none');


for j=1:n,
    
    if type(j)==1
        i=i+1;
        
        U(i)=Utemp(j);
        V(i)=Vtemp(j);
        VecU(i)=VecUtemp(j);
        VecV(i)=VecVtemp(j);
        
        try
            cellID(i)=strcat(num2str(filename),'_',map.cellID(j));
            pResponse{i,1}=map.Presponse{j};
            onset1(i,1)=map.onset1(j);
            OSI(i)=map.OSI(j);
            OSIBars(i)=map.OSIBars(j);
            DSI(i)=map.DSI(j);
            symmetry_ratio(i)=map.symmetry_ratio(j);
            sumCorrRep2(i)=map.sumCorrRep2(j);
            oriR2(i)=map.oriR2(j);
            poriSz(i)=map.poriSz(j);
            poriSz_bars(i)=map.poriSz_bars(j);
            poriSz_bars_ON(i)=map.poriSz_bars_ON(j);  
            poriSz_bars_OFF(i)=map.poriSz_bars_OFF(j);  
            original_rmse(i)=map.original_rmse(j);
            shuffle_rmse(i,1:100)=map.shuffle_rmse(j,:);
            theta(i,1:9)=map.theta(j,:);
            meanr(i,1:8)=map.meanr(j,:);
            meanr_all{i, 1}=map.meanr_all{j, 1};    
            ONorONOFForOFF(i)=map.ONorONOFForOFF(j);
            isVertical(i)=map.isVertical(j);
            isHorizontal(i)=map.isHorizontal(j);
            isCART(i)=map.isCART(j);
            isRBPMS(i)=map.isRBPMS(j);
            isRetro(i)=map.isRetro(j);
            if isRetro(i)
                if iscell(map.injectionSite)
                    injectionSiteTmp=[];
                    if cell2mat(strfind(map.injectionSite,'MTNd')); injectionSiteTmp=[injectionSiteTmp,1]; end
                    if cell2mat(strfind(map.injectionSite,'SF')); injectionSiteTmp=[injectionSiteTmp,2]; end
                    if cell2mat(strfind(map.injectionSite,'NOT')); injectionSiteTmp=[injectionSiteTmp,3]; end
                    if cell2mat(strfind(map.injectionSite,'IF')); injectionSiteTmp=[injectionSiteTmp,4]; end
                    injectionSite{i}=injectionSiteTmp;
                elseif ~iscell(map.injectionSite)
                    injectionSiteTmp=[];
                    if strfind(map.injectionSite,'MTNd'); injectionSiteTmp=[injectionSiteTmp,1]; end
                    if strfind(map.injectionSite,'SF'); injectionSiteTmp=[injectionSiteTmp,2]; end
                    if strfind(map.injectionSite,'NOT'); injectionSiteTmp=[injectionSiteTmp,3]; end
                    if strfind(map.injectionSite,'IF'); injectionSiteTmp=[injectionSiteTmp,4]; end
                    injectionSite{i}=injectionSiteTmp;
                end
            elseif isRetro(i)==0
                injectionSite{i}=[];
            end
        end
        
        %Interpolation of data onto grid
        
        if Sec(j)==1,
            
            
            %UV Data
            
            VecUcomp(i)=griddata(U1,V1,VecU1,U(i),V(i));
            VecVcomp(i)=griddata(U1,V1,VecV1,U(i),V(i));
            Ucomp(i)=U(i);
            Vcomp(i)=V(i);
            
            
            %STH data
            S(i)=griddata(U1,V1,S1,U(i),V(i));
            TH(i)=griddata(U1,V1,TH1,U(i),V(i));
            
            %Gen STH data
            Xprime(i)=griddata(U1,V1,Xprime1,U(i),V(i));
            Yprime(i)=griddata(U1,V1,Yprime1,U(i),V(i));
            
            %XYZ data
            X(i)=griddata(U1,V1,X1,U(i),V(i));
            Y(i)=griddata(U1,V1,Y1,U(i),V(i));
            Z(i)=griddata(U1,V1,Z1,U(i),V(i));
            
            %Conversion to STH vector field
            RHO(i)=griddata(S1,TH1,RHO1,S(i),TH(i));
            F(i)=griddata(S1,TH1,F1,S(i),TH(i));
            RHOS(i)=griddata(S1,TH1,RHOS1,S(i),TH(i));
            RHOTH(i)=griddata(S1,TH1,RHOTH1,S(i),TH(i));
            FS(i)=griddata(S1,TH1,FS1,S(i),TH(i));
            FTH(i)=griddata(S1,TH1,FTH1,S(i),TH(i));
            
            [g1,g2] = InverseVec(VecU(i),VecV(i),RHO(i),F(i),RHOS(i),FS(i),RHOTH(i),FTH(i));
  
                VecS(i)=g1;
                VecTH(i)=g2;
                
                VecScomp(i)=griddata(S1,TH1,VecS1,S(i),TH(i));
                VecTHcomp(i)=griddata(S1,TH1,VecTH1,S(i),TH(i));
                
                %Conversion to Generalized STH vector field
                
                VecXprime(i)=cos(TH(i)).*VecS(i)-S(i).*sin(TH(i)).*VecTH(i);
                VecYprime(i)=sin(TH(i)).*VecS(i)+S(i).*cos(TH(i)).*VecTH(i);
                
                VecXprimecomp(i)=griddata(Xprime1,Yprime1,VecXprime1,Xprime(i),Yprime(i));
                VecYprimecomp(i)=griddata(Xprime1,Yprime1,VecYprime1,Xprime(i),Yprime(i));
                
                %Conversion to Vectors on the sphere.
                VecX(i)=VecS(i).*cos(S(i)/R).*cos(TH(i))-R*VecTH(i).*sin(S(i)/R).*sin(TH(i));
                VecY(i)=VecS(i).*cos(S(i)/R).*sin(TH(i))+R*VecTH(i).*sin(S(i)/R).*cos(TH(i));
                VecZ(i)=VecS(i).*sin(S(i)/R);
                
                VecXcomp(i)=griddata(S1,TH1,VecX1,S(i),TH(i));
                VecYcomp(i)=griddata(S1,TH1,VecY1,S(i),TH(i));
                VecZcomp(i)=griddata(S1,TH1,VecZ1,S(i),TH(i));
                
        elseif Sec(j)==2,
            
            %UV Data
            
            VecUcomp(i)=griddata(U2,V2,VecU2,U(i),V(i));
            VecVcomp(i)=griddata(U2,V2,VecV2,U(i),V(i));
            Ucomp(i)=U(i);
            Vcomp(i)=V(i);
            
            %STH data
            S(i)=griddata(U2,V2,S2,U(i),V(i));
            TH(i)=griddata(U2,V2,TH2,U(i),V(i));
            
            %Gen STH data
            Xprime(i)=griddata(U2,V2,Xprime2,U(i),V(i));
            Yprime(i)=griddata(U2,V2,Yprime2,U(i),V(i));
            
            %XYZ data
            X(i)=griddata(U2,V2,X2,U(i),V(i));
            Y(i)=griddata(U2,V2,Y2,U(i),V(i));
            Z(i)=griddata(U2,V2,Z2,U(i),V(i));
            
            %Conversion to STH vector field
            RHO(i)=griddata(S2,TH2,RHO2,S(i),TH(i));
            F(i)=griddata(S2,TH2,F2,S(i),TH(i));
            RHOS(i)=griddata(S2,TH2,RHOS2,S(i),TH(i));
            RHOTH(i)=griddata(S2,TH2,RHOTH2,S(i),TH(i));
            FS(i)=griddata(S2,TH2,FS2,S(i),TH(i));
            FTH(i)=griddata(S2,TH2,FTH2,S(i),TH(i));
            
            [g1,g2] = InverseVec(VecU(i),VecV(i),RHO(i),F(i),RHOS(i),FS(i),RHOTH(i),FTH(i));

                VecS(i)=g1;
                VecTH(i)=g2;
                
                VecScomp(i)=griddata(S2,TH2,VecS2,S(i),TH(i));
                VecTHcomp(i)=griddata(S2,TH2,VecTH2,S(i),TH(i));
                
                %Conversion to Generalized STH vector field
                
                VecXprime(i)=cos(TH(i)).*VecS(i)-S(i).*sin(TH(i)).*VecTH(i);
                VecYprime(i)=sin(TH(i)).*VecS(i)+S(i).*cos(TH(i)).*VecTH(i);
                
                VecXprimecomp(i)=griddata(Xprime2,Yprime2,VecXprime2,Xprime(i),Yprime(i));
                VecYprimecomp(i)=griddata(Xprime2,Yprime2,VecYprime2,Xprime(i),Yprime(i));
                
                %Conversion to Vectors on the sphere.
                VecX(i)=VecS(i).*cos(S(i)/R).*cos(TH(i))-R*VecTH(i).*sin(S(i)/R).*sin(TH(i));
                VecY(i)=VecS(i).*cos(S(i)/R).*sin(TH(i))+R*VecTH(i).*sin(S(i)/R).*cos(TH(i));
                VecZ(i)=VecS(i).*sin(S(i)/R);
                
                VecXcomp(i)=griddata(S2,TH2,VecX2,S(i),TH(i));
                VecYcomp(i)=griddata(S2,TH2,VecY2,S(i),TH(i));
                VecZcomp(i)=griddata(S2,TH2,VecZ2,S(i),TH(i));
                
        elseif Sec(j)==3,
            
            %UV Data
            VecUcomp(i)=griddata(U3,V3,VecU3,U(i),V(i));
            VecVcomp(i)=griddata(U3,V3,VecV3,U(i),V(i));
            Ucomp(i)=U(i);
            Vcomp(i)=V(i);
            
            %STH data
            S(i)=griddata(U3,V3,S3,U(i),V(i));
            TH(i)=griddata(U3,V3,TH3,U(i),V(i));
            
            %Gen STH data
            Xprime(i)=griddata(U3,V3,Xprime3,U(i),V(i));
            Yprime(i)=griddata(U3,V3,Yprime3,U(i),V(i));
            
            %XYZ data
            X(i)=griddata(U3,V3,X3,U(i),V(i));
            Y(i)=griddata(U3,V3,Y3,U(i),V(i));
            Z(i)=griddata(U3,V3,Z3,U(i),V(i));
            
            %Conversion to STH vector field
            RHO(i)=griddata(S3,TH3,RHO3,S(i),TH(i));
            F(i)=griddata(S3,TH3,F3,S(i),TH(i));
            RHOS(i)=griddata(S3,TH3,RHOS3,S(i),TH(i));
            RHOTH(i)=griddata(S3,TH3,RHOTH3,S(i),TH(i));
            FS(i)=griddata(S3,TH3,FS3,S(i),TH(i));
            FTH(i)=griddata(S3,TH3,FTH3,S(i),TH(i));
            
            [g1,g2] = InverseVec(VecU(i),VecV(i),RHO(i),F(i),RHOS(i),FS(i),RHOTH(i),FTH(i));
            
                VecS(i)=g1;
                VecTH(i)=g2;
                
                VecScomp(i)=griddata(S3,TH3,VecS3,S(i),TH(i));
                VecTHcomp(i)=griddata(S3,TH3,VecTH3,S(i),TH(i));
                
                %Conversion to Generalized STH vector field
                
                VecXprime(i)=cos(TH(i)).*VecS(i)-S(i).*sin(TH(i)).*VecTH(i);
                VecYprime(i)=sin(TH(i)).*VecS(i)+S(i).*cos(TH(i)).*VecTH(i);
                
                VecXprimecomp(i)=griddata(Xprime3,Yprime3,VecXprime3,Xprime(i),Yprime(i));
                VecYprimecomp(i)=griddata(Xprime3,Yprime3,VecYprime3,Xprime(i),Yprime(i));
                
                %Conversion to Vectors on the sphere.
                VecX(i)=VecS(i).*cos(S(i)/R).*cos(TH(i))-R*VecTH(i).*sin(S(i)/R).*sin(TH(i));
                VecY(i)=VecS(i).*cos(S(i)/R).*sin(TH(i))+R*VecTH(i).*sin(S(i)/R).*cos(TH(i));
                VecZ(i)=VecS(i).*sin(S(i)/R);
                
                VecXcomp(i)=griddata(S3,TH3,VecX3,S(i),TH(i));
                VecYcomp(i)=griddata(S3,TH3,VecY3,S(i),TH(i));
                VecZcomp(i)=griddata(S3,TH3,VecZ3,S(i),TH(i));
                
        elseif Sec(j)==4
            
            %UV Data
            VecUcomp(i)=griddata(U4,V4,VecU4,U(i),V(i));
            VecVcomp(i)=griddata(U4,V4,VecV4,U(i),V(i));
            Ucomp(i)=U(i);
            Vcomp(i)=V(i);
            
            %STH data
            S(i)=griddata(U4,V4,S4,U(i),V(i));
            TH(i)=griddata(U4,V4,TH4,U(i),V(i));
            
            %Gen STH data
            Xprime(i)=griddata(U4,V4,Xprime4,U(i),V(i));
            Yprime(i)=griddata(U4,V4,Yprime4,U(i),V(i));
            
            %XYZ data
            X(i)=griddata(U4,V4,X4,U(i),V(i));
            Y(i)=griddata(U4,V4,Y4,U(i),V(i));
            Z(i)=griddata(U4,V4,Z4,U(i),V(i));
            
            %Conversion to STH vector field
            RHO(i)=griddata(S4,TH4,RHO4,S(i),TH(i));
            F(i)=griddata(S4,TH4,F4,S(i),TH(i));
            RHOS(i)=griddata(S4,TH4,RHOS4,S(i),TH(i));
            RHOTH(i)=griddata(S4,TH4,RHOTH4,S(i),TH(i));
            FS(i)=griddata(S4,TH4,FS4,S(i),TH(i));
            FTH(i)=griddata(S4,TH4,FTH4,S(i),TH(i));
            
            [g1,g2] = InverseVec(VecU(i),VecV(i),RHO(i),F(i),RHOS(i),FS(i),RHOTH(i),FTH(i));
                  
                VecS(i)=g1;
                VecTH(i)=g2;
                
                VecScomp(i)=griddata(S4,TH4,VecS4,S(i),TH(i));
                VecTHcomp(i)=griddata(S4,TH4,VecTH4,S(i),TH(i));
                
                %Conversion to Generalized STH vector field
                
                VecXprime(i)=cos(TH(i)).*VecS(i)-S(i).*sin(TH(i)).*VecTH(i);
                VecYprime(i)=sin(TH(i)).*VecS(i)+S(i).*cos(TH(i)).*VecTH(i);
                
                VecXprimecomp(i)=griddata(Xprime4,Yprime4,VecXprime4,Xprime(i),Yprime(i));
                VecYprimecomp(i)=griddata(Xprime4,Yprime4,VecYprime4,Xprime(i),Yprime(i));
                
                %Conversion to Vectors on the sphere.
                VecX(i)=VecS(i).*cos(S(i)/R).*cos(TH(i))-R*VecTH(i).*sin(S(i)/R).*sin(TH(i));
                VecY(i)=VecS(i).*cos(S(i)/R).*sin(TH(i))+R*VecTH(i).*sin(S(i)/R).*cos(TH(i));
                VecZ(i)=VecS(i).*sin(S(i)/R);
                
                VecXcomp(i)=griddata(S4,TH4,VecX4,S(i),TH(i));
                VecYcomp(i)=griddata(S4,TH4,VecY4,S(i),TH(i));
                VecZcomp(i)=griddata(S4,TH4,VecZ4,S(i),TH(i));
                
        end
    end
end


%% Normalizing Vector Fields
temp=(VecU.^2+VecV.^2).^(1/2);
VecU=VecU./temp;
VecV=VecV./temp;

temp=(VecUcomp.^2+VecVcomp.^2).^(1/2);
VecUcomp=VecUcomp./temp;
VecVcomp=VecVcomp./temp;


temp=(VecS.^2+VecTH.^2).^(1/2);
VecS=VecS./temp;
VecTH=VecTH./temp;

temp=(VecScomp.^2+VecTHcomp.^2).^(1/2);
VecScomp=VecScomp./temp;
VecTHcomp=VecTHcomp./temp;

ind=isnan(VecS);
S(ind)=[];
TH(ind)=[];
VecS(ind)=[];
VecTH(ind)=[];
VecScomp(ind)=[];
VecTHcomp(ind)=[];

U(ind)=[];
V(ind)=[];
VecU(ind)=[];
VecV(ind)=[];
VecUcomp(ind)=[];
VecVcomp(ind)=[];
Ucomp(ind)=[];
Vcomp(ind)=[];

temp=(VecXprime.^2+VecYprime.^2).^(1/2);
VecXprime=VecXprime./temp;
VecYprime=VecYprime./temp;

temp=(VecXprimecomp.^2+VecYprimecomp.^2).^(1/2);
VecXprimecomp=VecXprimecomp./temp;
VecYprimecomp=VecYprimecomp./temp;

ind=isnan(VecXprime);
Xprime(ind)=[];
Yprime(ind)=[];
VecXprime(ind)=[];
VecYprime(ind)=[];
VecXprimecomp(ind)=[];
VecYprimecomp(ind)=[];

temp=(VecX.^2+VecY.^2+VecZ.^2).^(1/2);
VecX=VecX./temp;
VecY=VecY./temp;
VecZ=VecZ./temp;

temp=(VecXcomp.^2+VecYcomp.^2+VecZcomp.^2).^(1/2);
VecXcomp=VecXcomp./temp;
VecYcomp=VecYcomp./temp;
VecZcomp=VecZcomp./temp;

% removing cells right on the singularity
ind=isnan(VecX);
ind2=find(ind==1);
X(ind)=[];
Y(ind)=[];
Z(ind)=[];
VecX(ind)=[];
VecY(ind)=[];
VecZ(ind)=[];
VecXcomp(ind)=[];
VecYcomp(ind)=[];
VecZcomp(ind)=[];

for g=numel(ind2):-1:1
cellID(ind2(g))=[];
injectionSite(ind2(g))=[];
pResponse(ind2(g))=[];
end
isCART(ind)=[];
isRBPMS(ind)=[];
isRetro(ind)=[];
isVertical(ind)=[];
isHorizontal(ind)=[];
ONorONOFForOFF(ind)=[];
OSI(ind)=[];
OSIBars(ind)=[];
DSI(ind)=[];
symmetry_ratio(ind)=[];
sumCorrRep2(ind)=[];
oriR2(ind)=[];
poriSz(ind)=[];
poriSz_bars(ind)=[];
poriSz_bars_ON(ind)=[];  
poriSz_bars_OFF(ind)=[];
original_rmse(ind)=[];
original_rmse=original_rmse';
shuffle_rmse(ind',:)=[];
meanr(ind',:)=[];
meanr_all(ind)=[];
theta(ind',:)=[];
onset1(ind)=[];

for i=1:size(meanr,1)
a=meanr(i,:);
[~,I]=max(a);
meanrs(i,1:8)=circshift(a,-(I-1));
end


%% Computing Angles Between Vectors

angleUV=acos(VecU.*VecUcomp+VecV.*VecVcomp);
angleSTH=acos(VecS.*VecScomp+VecTH.*VecTHcomp);
angleSTHGen=acos(VecXprime.*VecXprimecomp+VecYprime.*VecYprimecomp);
angleXYZ=acos(VecX.*VecXcomp+VecY.*VecYcomp+VecZ.*VecZcomp);

% figure;
% hold on;
% plot(angleUV,'o');
% plot(angleSTH,'o');
% plot(angleSTHGen,'o');
% plot(angleXYZ,'o');


UVvecComp.U=U;
UVvecComp.V=V;
UVvecComp.VecU=VecU;
UVvecComp.VecV=VecV;
UVvecComp.VecUcomp=VecUcomp;
UVvecComp.VecVcomp=VecVcomp;
UVvecComp.angleUV=angleUV;

STHvecComp.S=S;
STHvecComp.TH=TH;
STHvecComp.VecS=VecS;
STHvecComp.VecTH=VecTH;
STHvecComp.VecScomp=VecScomp;
STHvecComp.VecTHcomp=VecTHcomp;
STHvecComp.angleSTH=angleSTH;

STHgenvecComp.Xprime=Xprime;
STHgenvecComp.Yprime=Yprime;
STHgenvecComp.VecXprime=VecXprime;
STHgenvecComp.VecYprime=VecYprime;
STHgenvecComp.VecXprimecomp=VecXprimecomp;
STHgenvecComp.VecYprimecomp=VecYprimecomp;
STHgenvecComp.angleSTHGen=angleSTHGen;

XYZvecComp.X=X;
XYZvecComp.Y=Y;
XYZvecComp.Z=Z;
XYZvecComp.VecX=VecX;
XYZvecComp.VecY=VecY;
XYZvecComp.VecZ=VecZ;
XYZvecComp.VecXcomp=VecXcomp;
XYZvecComp.VecYcomp=VecYcomp;
XYZvecComp.VecZcomp=VecZcomp;
XYZvecComp.angleXYZ=angleXYZ;



%Plotting vector fields in UV coordinates
figure;
hold on;
surf(U1,V1,zeros(size(U1,1),size(U1,1)), 'EdgeColor', 'none');
surf(U2,V2,zeros(size(U2,1),size(U2,1)), 'EdgeColor', 'none');
surf(U3,V3,zeros(size(U3,1),size(U3,1)), 'EdgeColor', 'none');
surf(U4,V4,zeros(size(U4,1),size(U4,1)), 'EdgeColor', 'none');
scale=0.4;
%quiver(Ucomp,Vcomp,VecUcomp,VecVcomp,scale);
quiver(U,V,VecU,VecV,scale);
drawnow expose
hold off;


%Plotting vector fields in STH coordinates
% figure;
% hold on;
% scale=0.2;
% quiver(S,TH/(2*pi),VecScomp,VecTHcomp,scale);
% quiver(S,TH/(2*pi),VecS,VecTH,scale);
% drawnow expose
% hold off;

%Plotting vector fields in generalized (S,TH) coordinates
% figure;
% hold on;
% scale=0.2;
% quiver(Xprime,Yprime,VecXprimecomp,VecYprimecomp,scale);
% quiver(Xprime,Yprime,VecXprime,VecYprime,scale);
% drawnow expose
% hold off;

%Plotting vector fields in generalized (X,Y,Z) coordinates
% figure;
% hold on;
% scale=0.2;
% quiver3(X,Y,Z,VecXcomp,VecYcomp,VecZcomp,scale);
% quiver3(X,Y,Z,VecX,VecY,VecZ,scale);
% drawnow expose
% hold off;


MaxIndex=i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Mapping to "standard" retina.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

n=size(StSec1Data,1);

StRHO1=StSec1Data(1:n,1:n);
StF1=StSec1Data(1:n,n+1:2*n);
StS1=StSec1Data(1:n,2*n+1:3*n);
StTH1=StSec1Data(1:n,3*n+1:4*n);

StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);

%Calcution of derivatives
ds=StS1(1,2)-StS1(1,1);
DS=DiffX(n)/ds;
StRHO1S=(DS*StRHO1')';
StF1S=(DS*StF1')';

dth=StTH1(2,1)-StTH1(1,1);
DTH=DiffX(n)/dth;
StRHO1TH=(DTH*StRHO1);
StF1TH=(DTH*StF1);

StRHO2=StSec2Data(1:n,1:n);
StF2=StSec2Data(1:n,n+1:2*n);
StS2=StSec2Data(1:n,2*n+1:3*n);
StTH2=StSec2Data(1:n,3*n+1:4*n);

StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);

%Calcution of derivatives
ds=StS2(1,2)-StS2(1,1);
DS=DiffX(n)/ds;
StRHO2S=(DS*StRHO2')';
StF2S=(DS*StF2')';

dth=StTH2(2,1)-StTH2(1,1);
DTH=DiffX(n)/dth;
StRHO2TH=(DTH*StRHO2);
StF2TH=(DTH*StF2);

StRHO3=StSec3Data(1:n,1:n);
StF3=StSec3Data(1:n,n+1:2*n);
StS3=StSec3Data(1:n,2*n+1:3*n);
StTH3=StSec3Data(1:n,3*n+1:4*n);

StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);

%Calcution of derivatives
ds=StS3(1,2)-StS3(1,1);
DS=DiffX(n)/ds;
StRHO3S=(DS*StRHO3')';
StF3S=(DS*StF3')';

dth=StTH3(2,1)-StTH3(1,1);
DTH=DiffX(n)/dth;
StRHO3TH=(DTH*StRHO3);
StF3TH=(DTH*StF3);

StRHO4=StSec4Data(1:n,1:n);
StF4=StSec4Data(1:n,n+1:2*n);
StS4=StSec4Data(1:n,2*n+1:3*n);
StTH4=StSec4Data(1:n,3*n+1:4*n);

StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);

%Calcution of derivatives
ds=StS4(1,2)-StS4(1,1);
DS=DiffX(n)/ds;
StRHO4S=(DS*StRHO4')';
StF4S=(DS*StF4')';

dth=StTH4(2,1)-StTH4(1,1);
DTH=DiffX(n)/dth;
StRHO4TH=(DTH*StRHO4);
StF4TH=(DTH*StF4);

for i=1:size(TH,2)
    angle=mod(TH(i),2*pi);
    if angle<pi/2
        
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i)-StS1).^2+(angle-StTH1).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i)=StRHO1(a(i),b(i))*cos(StF1(a(i),b(i)));
        VSt(i)=StRHO1(a(i),b(i))*sin(StF1(a(i),b(i)));
        
        VecUSt(i)=VecS(i).*(StRHO1S(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*sin(StF1(a(i),b(i))))+VecTH(i).*(StRHO1TH(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*sin(StF1(a(i),b(i))));
        VecVSt(i)=VecS(i).*(StRHO1S(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*cos(StF1(a(i),b(i))))+VecTH(i).*(StRHO1TH(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*cos(StF1(a(i),b(i))));
        
        VecUStcomp(i)=VecScomp(i).*(StRHO1S(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*sin(StF1(a(i),b(i))))+VecTHcomp(i).*(StRHO1TH(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*sin(StF1(a(i),b(i))));
        VecVStcomp(i)=VecScomp(i).*(StRHO1S(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*cos(StF1(a(i),b(i))))+VecTHcomp(i).*(StRHO1TH(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*cos(StF1(a(i),b(i))));
        
    elseif angle<pi
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i)-StS2).^2+(angle-StTH2).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i)=StRHO2(a(i),b(i))*cos(StF2(a(i),b(i)));
        VSt(i)=StRHO2(a(i),b(i))*sin(StF2(a(i),b(i)));
        
        VecUSt(i)=VecS(i).*(StRHO2S(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*sin(StF2(a(i),b(i))))+VecTH(i).*(StRHO2TH(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*sin(StF2(a(i),b(i))));
        VecVSt(i)=VecS(i).*(StRHO2S(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*cos(StF2(a(i),b(i))))+VecTH(i).*(StRHO2TH(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*cos(StF2(a(i),b(i))));
        
        VecUStcomp(i)=VecScomp(i).*(StRHO2S(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*sin(StF2(a(i),b(i))))+VecTHcomp(i).*(StRHO2TH(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*sin(StF2(a(i),b(i))));
        VecVStcomp(i)=VecScomp(i).*(StRHO2S(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*cos(StF2(a(i),b(i))))+VecTHcomp(i).*(StRHO2TH(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*cos(StF2(a(i),b(i))));
        
    elseif angle<3*pi/2
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i)-StS3).^2+(angle-StTH3).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i)=StRHO3(a(i),b(i))*cos(StF3(a(i),b(i)));
        VSt(i)=StRHO3(a(i),b(i))*sin(StF3(a(i),b(i)));
        
        VecUSt(i)=VecS(i).*(StRHO3S(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*sin(StF3(a(i),b(i))))+VecTH(i).*(StRHO3TH(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*sin(StF3(a(i),b(i))));
        VecVSt(i)=VecS(i).*(StRHO3S(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*cos(StF3(a(i),b(i))))+VecTH(i).*(StRHO3TH(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*cos(StF3(a(i),b(i))));
        
        VecUStcomp(i)=VecScomp(i).*(StRHO3S(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*sin(StF3(a(i),b(i))))+VecTHcomp(i).*(StRHO3TH(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*sin(StF3(a(i),b(i))));
        VecVStcomp(i)=VecScomp(i).*(StRHO3S(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*cos(StF3(a(i),b(i))))+VecTHcomp(i).*(StRHO3TH(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*cos(StF3(a(i),b(i))));
        
    else
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i)-StS4).^2+(angle-StTH4).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i)=StRHO4(a(i),b(i))*cos(StF4(a(i),b(i)));
        VSt(i)=StRHO4(a(i),b(i))*sin(StF4(a(i),b(i)));
        
        VecUSt(i)=VecS(i).*(StRHO4S(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*sin(StF4(a(i),b(i))))+VecTH(i).*(StRHO4TH(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*sin(StF4(a(i),b(i))));
        VecVSt(i)=VecS(i).*(StRHO4S(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*cos(StF4(a(i),b(i))))+VecTH(i).*(StRHO4TH(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*cos(StF4(a(i),b(i))));
        
        VecUStcomp(i)=VecScomp(i).*(StRHO4S(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*sin(StF4(a(i),b(i))))+VecTHcomp(i).*(StRHO4TH(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*sin(StF4(a(i),b(i))));
        VecVStcomp(i)=VecScomp(i).*(StRHO4S(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*cos(StF4(a(i),b(i))))+VecTHcomp(i).*(StRHO4TH(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*cos(StF4(a(i),b(i))));
        
    end
end

%figure;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

temp=(VecUStcomp.^2+VecVStcomp.^2).^(1/2);
VecUStcomp=VecUStcomp./temp;
VecVStcomp=VecVStcomp./temp;


%Plotting vector fields in Standard UV coordinates
figure;
hold on;
surf(StU1,StV1,zeros(size(StU1,1),size(StU1,1)), 'EdgeColor', 'none');
surf(StU2,StV2,zeros(size(StU2,1),size(StU2,1)), 'EdgeColor', 'none');
surf(StU3,StV3,zeros(size(StU3,1),size(StU3,1)), 'EdgeColor', 'none');
surf(StU4,StV4,zeros(size(StU4,1),size(StU4,1)), 'EdgeColor', 'none');
scale=0.4;

quiver(USt,VSt,VecUStcomp,VecVStcomp,scale);
quiver(USt,VSt,VecUSt,VecVSt,scale);
drawnow expose
hold off;

angleUVSt=acos(VecUSt.*VecUStcomp+VecVSt.*VecVStcomp);


UVStvecComp.USt=USt;
UVStvecComp.VSt=VSt;
UVStvecComp.VecUSt=VecUSt;
UVStvecComp.VecVSt=VecVSt;
UVStvecComp.VecUStcomp=VecUStcomp;
UVStvecComp.VecVStcomp=VecVStcomp;
UVStvecComp.angleUVSt=angleUVSt;

try
    save(['vecCompAlphaCorrBars_OFF_',num2str(discPoints),'_','DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_',fieldType,'_',filename,'.mat'],'UVvecComp','STHvecComp','STHgenvecComp',...
        'XYZvecComp','UVStvecComp','ONorONOFForOFF','isRetro','injectionSite','isCART','isRBPMS','isVertical','isHorizontal','cellID','DSI','OSI','symmetry_ratio','sumCorrRep2','oriR2','poriSz','poriSz_bars','poriSz_bars_ON','poriSz_bars_OFF','original_rmse','shuffle_rmse','meanrs','meanr','meanr_all','OSIBars','theta','pResponse','onset1');
catch
    save(['vecCompAlphaCorrBars_OFF_',num2str(discPoints),'_','DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_',fieldType,'_',filename,'.mat'],'UVvecComp','STHvecComp','STHgenvecComp','XYZvecComp','UVStvecComp');
end
% dlmwrite('angleStUVdif',angleStUV);



end


