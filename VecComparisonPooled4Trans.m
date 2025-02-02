function [allPossibleFields]= VecComparisonPooled4Trans(data,a,maxVec,alpha,beta,polarity,fieldID,VizAngle,allPossibleFields,k,h,interpType,histType)
                                                      
% histType can be either 1 (predicted vector for ASC), 2 (predicted vector
% for PSC), 3 (predicted vector for LSC), 4 2 (predicted vector for all
% three canals), or 5 (real preferred direction of cells).

% cellType can be either 'ONDS', 'ONOFFDS', 'OFFDS', or 'retroONDS'.
% fieldID can be either 'ASC', 'PSC', 'LSC', 'cardinalUp', etc.
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

%% Test Parameters Comment out when not testing.
% VizAngle=180;
% normalASC=[0.592,0.764,0.247];     %ASC
% [alpha,beta]=AB(normalASC);
% polarity=1;
%% Physical Parameters.
% This data is imported from data exported from the mapping script. These
% values correspond to measured values on the flattened retina.

phys=dlmread('phys');

R=phys(1);

M=phys(2); %This should be normalized to 1.

VizAngle=deg2rad(VizAngle);

switch interpType
    case 1          % perform interpolation
        n=40;       % Number of points in discretization of vector field on the sphere.
    case 2          % find closest grid point
        n=100;
end

%% Creating Vector field in XYZ coordinates.
% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(0,2*pi,n);


[S,TH]=meshgrid(s,th);
%Calculation of the vector field of rotations in the coordinates s and th.
%The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation.

%% Creating Vector field in Visual Space coordinates.
SV=-VizAngle*R/(2*M)*S+pi*R;
THV=TH+pi;

t=zeros(n,n);
phi=zeros(n,n);


% normalASC=[0.592,0.764,0.247];     %ASC
% [alpha,beta]=AB(normalASC);



%This loop constructs the components of the vector field of rotations
%in the the s, th coordinates. The function InverseSTH determines the
%coordinates of (t,phi) given coordinates s and th and is necessary to
%calculate the components of the vector field.

for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,SV(i,j),THV(i,j));
        
        vecSV(i,j)=sin(SV(i,j)/R)*(cos(beta)*sin(t(i,j)/R)+...
            cos(t(i,j)/R)*cos(phi(i,j))*sin(beta));
        vecSV(i,j)=vecSV(i,j)+cos(SV(i,j)/R)*(-cos(alpha-THV(i,j))*sin(t(i,j)/R)*sin(beta)+...
            cos(t(i,j)/R)*(cos(beta)*cos(alpha-THV(i,j))*cos(phi(i,j))-...
            sin(alpha-THV(i,j))*sin(phi(i,j))));
        
        vecTHV(i,j)=-sin(t(i,j)/R)*sin(beta)*sin(alpha-THV(i,j))+cos(t(i,j)/R)*...
            (cos(beta)*cos(phi(i,j))*sin(alpha-THV(i,j))+cos(alpha-THV(i,j))*sin(phi(i,j)));
        vecTHV(i,j)=1/R*csc(SV(i,j)/R)*vecTHV(i,j);
        
    end
end


temp=(vecSV.^2+vecTHV.^2).^(1/2);

vecSV=vecSV./temp;
vecTHV=vecTHV./temp;

% figure;
% hold on;
% quiver(SV,THV,vecSV,vecTHV,.2);
% plot(R*beta,alpha,'o');
% hold off;

% figure;
XV=R*sin(SV/R).*cos(THV);
YV=R*sin(SV/R).*sin(THV);
ZV=-R*cos(SV/R);

vecXV=vecSV.*cos(SV/R).*cos(THV)-R*vecTHV.*sin(SV/R).*sin(THV);
vecYV=vecSV.*cos(SV/R).*sin(THV)+R*vecTHV.*sin(SV/R).*cos(THV);
vecZV=vecSV.*sin(SV/R);

%quiver3(XV,YV,ZV,vecXV,vecYV,vecZV);

temp=(vecSV.^2+vecTHV.^2).^(1/2);
vecSV=vecSV./temp;
vecTHV=vecTHV./temp;

vecS=-2*M/(VizAngle*R)*vecSV;

vecTH=vecTHV.*polarity;

% figure;
% quiver(S,TH,vecS,vecTH,.2);
% hold on;
% plot(R*beta,alpha,'o');
% hold off;

discPoints=40;


vecX=vecS.*cos(S/R).*cos(TH)-R*vecTH.*sin(S/R).*sin(TH);
vecY=vecS.*cos(S/R).*sin(TH)+R*vecTH.*sin(S/R).*cos(TH);
vecZ=vecS.*sin(S/R);

% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

XM=R.*cos(TH).*sin(S./R);
YM=R.*sin(TH).*sin(S./R);
ZM=R-R.*cos(S./R);

vecXM=veca;
vecYM=vecb;
vecZM=vecc;

% figure;
%quiver3(XM,YM,ZM,vecXM,vecYM,vecZM);


%% Input measured data mapped to XYZ coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   We now input measured Data from the retina in the flattened
%   coordinates.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% import XYZ data.

switch histType
    case 5
        %%% observed vectors are set to the preferred direction of cells
        %%% (real)
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=data.XYZvecComp.VecX';
        VecY=data.XYZvecComp.VecY';
        VecZ=data.XYZvecComp.VecZ';
        
%         load('maxVec_Calcium imaging data_translation_25percentAll4ONOFF_40_DS_created_Jul_01_2015_19_16.mat');
%         VecX=maxVec.X';
%         VecY=maxVec.Y';
%         VecZ=maxVec.Z';
        
  
    case 6
        %%% observed vectors are set to the preferred direction of cells
        %%% from combinedMao (fake)
        X=data.X';
        Y=data.Y';
        Z=data.Z';
        VecX=data.VecX';
        VecY=data.VecY';
        VecZ=data.VecZ';
        
        
        case 7
        %%% observed vectors are set to the preferred direction of cells
        %%% from a single-channel map (fake)
        X=data.X';
        Y=data.Y';
        Z=data.Z';
        VecX=data.VecX';
        VecY=data.VecY';
        VecZ=data.VecZ';
        
        
        
    case 1
        %%% observed vectors are set to the vector field around the ASC canal
        
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=data.XYZvecComp.VecXcompASC';
        VecY=data.XYZvecComp.VecYcompASC';
        VecZ=data.XYZvecComp.VecZcompASC';
        
    case 2
        %%% observed vectors are set to the vector field around the PSC
        %%% canal (inverted)
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=data.XYZvecComp.VecXcompPSC';
        VecY=data.XYZvecComp.VecYcompPSC';
        VecZ=data.XYZvecComp.VecZcompPSC';
        
    case 3
        %%% observed vectors are set to the vector field around the LSC canal
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=data.XYZvecComp.VecXcompLSC';
        VecY=data.XYZvecComp.VecYcompLSC';
        VecZ=data.XYZvecComp.VecZcompLSC';
        
    case -1
        %%% observed vectors are set to the vector field around the ASC canal
        
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=-data.XYZvecComp.VecXcompASC';
        VecY=-data.XYZvecComp.VecYcompASC';
        VecZ=-data.XYZvecComp.VecZcompASC';
        
    case -2
        %%% observed vectors are set to the vector field around the PSC canal
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=-data.XYZvecComp.VecXcompPSC';
        VecY=-data.XYZvecComp.VecYcompPSC';
        VecZ=-data.XYZvecComp.VecZcompPSC';
        
    case -3
        %%% observed vectors are set to the vector field around the LSC canal
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=-data.XYZvecComp.VecXcompLSC';
        VecY=-data.XYZvecComp.VecYcompLSC';
        VecZ=-data.XYZvecComp.VecZcompLSC';
        
        
    case 4
        %%% 50% of observed vectors are set to the vector field around the LSC and 50% to the INVERSE vector field around the LSC
        %         f=randi(2,[size(data.XYZvecComp.VecXcompLSC,1),1]);
        %         X=[data.XYZvecComp.X(f==1);data.XYZvecComp.X(f==2)];
        %         Y=[data.XYZvecComp.Y(f==1);data.XYZvecComp.Y(f==2)];
        %         Z=[data.XYZvecComp.Z(f==1);data.XYZvecComp.Z(f==2)];
        %         X=X';
        %         Y=Y';
        %         Z=Z';
        %
        %         VecX=[data.XYZvecComp.VecXcompLSC(f==1);-data.XYZvecComp.VecXcompLSC(f==2)];
        %         VecY=[data.XYZvecComp.VecYcompLSC(f==1);-data.XYZvecComp.VecYcompLSC(f==2)];
        %         VecZ=[data.XYZvecComp.VecZcompLSC(f==1);-data.XYZvecComp.VecZcompLSC(f==2)];
        %         VecX=VecX';
        %         VecY=VecY';
        %         VecZ=VecZ';
        
        %%% 25% ASC, 25% PSC, 25% LSC, and 25% inversePSC.
        f=randi(4,[size(data.XYZvecComp.VecXcompLSC,1),1]);
        X=[data.XYZvecComp.X(f==1);data.XYZvecComp.X(f==2);data.XYZvecComp.X(f==3);data.XYZvecComp.X(f==4)];
        Y=[data.XYZvecComp.Y(f==1);data.XYZvecComp.Y(f==2);data.XYZvecComp.Y(f==3);data.XYZvecComp.Y(f==4)];
        Z=[data.XYZvecComp.Z(f==1);data.XYZvecComp.Z(f==2);data.XYZvecComp.Z(f==3);data.XYZvecComp.Z(f==4)];
        X=X';
        Y=Y';
        Z=Z';
        
        VecX=[data.XYZvecComp.VecXcompASC(f==1);data.XYZvecComp.VecXcompPSC(f==2);data.XYZvecComp.VecXcompLSC(f==3);-data.XYZvecComp.VecXcompPSC(f==4)];
        VecY=[data.XYZvecComp.VecYcompASC(f==1);data.XYZvecComp.VecYcompPSC(f==2);data.XYZvecComp.VecYcompLSC(f==3);-data.XYZvecComp.VecYcompPSC(f==4)];
        VecZ=[data.XYZvecComp.VecZcompASC(f==1);data.XYZvecComp.VecZcompPSC(f==2);data.XYZvecComp.VecZcompLSC(f==3);-data.XYZvecComp.VecZcompPSC(f==4)];
        VecX=VecX';
        VecY=VecY';
        VecZ=VecZ';
        
        
        %         f=randi(3,[size(data.XYZvecComp.VecXcompLSC,1),1]);
        %         X=[data.XYZvecComp.X(f==1);data.XYZvecComp.X(f==2);data.XYZvecComp.X(f==3)];
        %         Y=[data.XYZvecComp.Y(f==1);data.XYZvecComp.Y(f==2);data.XYZvecComp.Y(f==3)];
        %         Z=[data.XYZvecComp.Z(f==1);data.XYZvecComp.Z(f==2);data.XYZvecComp.Z(f==3)];
        %         X=X';
        %         Y=Y';
        %         Z=Z';
        %
        %         VecX=[data.XYZvecComp.VecXcompASC(f==1);data.XYZvecComp.VecXcompPSC(f==2);data.XYZvecComp.VecXcompLSC(f==3)];
        %         VecY=[data.XYZvecComp.VecYcompASC(f==1);data.XYZvecComp.VecYcompPSC(f==2);data.XYZvecComp.VecYcompLSC(f==3)];
        %         VecZ=[data.XYZvecComp.VecZcompASC(f==1);data.XYZvecComp.VecZcompPSC(f==2);data.XYZvecComp.VecZcompLSC(f==3)];
        %         VecX=VecX';
        %         VecY=VecY';
        %         VecZ=VecZ';
        
    case 10
        
        X=data.XYZvecComp.X;
        Y=data.XYZvecComp.Y;
        Z=data.XYZvecComp.Z;
        X=X';
        Y=Y';
        Z=Z';
        
        % randomize vecX, vecY, and vecZ
        %          fx=-1+(1+1)*rand(size(data.XYZvecComp.VecXcompLSC,1),1);
        %          fy=-1+(1+1)*rand(size(data.XYZvecComp.VecXcompLSC,1),1);
        %          fz=-1+(1+1)*rand(size(data.XYZvecComp.VecXcompLSC,1),1);
        %
        %          VecX=fx;
        %          VecY=fy;
        %          VecZ=fz;
        %          VecX=VecX';
        %          VecY=VecY';
        %          VecZ=VecZ';
        
        % randomize angle a
        %         a=deg2rad(360*rand(size(data.XYZvecComp.VecXcompLSC,1),1));
        %         a=a';
        
        VecX=(1./(sqrt(X.^2+Y.^2))).*((((R-Z).*sin(a).*X)./R)-cos(a).*Y);
        VecY=(1./(sqrt(X.^2+Y.^2))).*((((R-Z).*sin(a).*Y)./R)+cos(a).*X);
        VecZ=(sqrt(2.*R.*Z-Z.^2).*sin(a))./R;
        
        
        
        %Normalizing Vector Fields
        
        temp1=(VecX.^2+VecY.^2+VecZ.^2).^(1/2);
        VecX=VecX./temp1;
        VecY=VecY./temp1;
        VecZ=VecZ./temp1;
        
    case 20
        %%% observed vectors are set to the preferred direction of cells
        %%% (real)
        X=data.XYZvecComp.X';
        Y=data.XYZvecComp.Y';
        Z=data.XYZvecComp.Z';
        VecX=real(maxVec.X');
        VecY=real(maxVec.Y');
        VecZ=real(maxVec.Z');
        
        
end
%% Interpolating vector field at data points

switch interpType
    
    case 1          % perform interpolation
        XMshaped=reshape(XM,n^2,1);
        YMshaped=reshape(YM,n^2,1);
        ZMshaped=reshape(ZM,n^2,1);
        VecXMshaped=reshape(vecXM,n^2,1);
        VecYMshaped=reshape(vecYM,n^2,1);
        VecZMshaped=reshape(vecZM,n^2,1);
        
      
        VecXMInterp=scatteredInterpolant(XMshaped,YMshaped,ZMshaped,VecXMshaped);
        VecYMInterp=scatteredInterpolant(XMshaped,YMshaped,ZMshaped,VecYMshaped);
        VecZMInterp=scatteredInterpolant(XMshaped,YMshaped,ZMshaped,VecZMshaped);
        
        % figure;
        % quiver3(XM,YM,ZM,vecXM,vecYM,vecZM);
        
        Npooled=size(X,2);
        for i=1:Npooled
            %XYZ data
            VecXcomp(i)=VecXMInterp(X(i),Y(i),Z(i));
            VecYcomp(i)=VecYMInterp(X(i),Y(i),Z(i));
            VecZcomp(i)=VecZMInterp(X(i),Y(i),Z(i));
        end
        
        % figure;
        % quiver3(X,Y,Z,VecX,VecY,VecZ);
        
    case 2      % find closest grid point
        
        Npooled=size(X,2);
        for i=1:Npooled
            temp=(X(i)-XM).^2+(Y(i)-YM).^2+(Z(i)-ZM).^2;
            [ind1,ind2] =find(temp==min(min(temp)));
            a(i)=ind1(1);
            b(i)=ind2(1);
            VecXcomp(i)=vecXM(a(i),b(i));
            VecYcomp(i)=vecYM(a(i),b(i));
            VecZcomp(i)=vecZM(a(i),b(i));
        end
        
    case 3
        
        VecXcomp=griddata(XM,YM,ZM,VecXM,X,Y,Z);
        VecYcomp=griddata(XM,YM,ZM,VecYM,X,Y,Z);
        VecZcomp=griddata(XM,YM,ZM,VecZM,X,Y,Z);
end

%% Normalizing Vector Fields

temp=(VecXcomp.^2+VecYcomp.^2+VecZcomp.^2).^(1/2);
VecXcomp=VecXcomp./temp;
VecYcomp=VecYcomp./temp;
VecZcomp=VecZcomp./temp;

% figure;
% quiver3(X,Y,Z,VecX,VecY,VecZ,'r');
% hold on
% quiver3(X,Y,Z,VecXcomp,VecYcomp,VecZcomp,'b');

%% calculate meanOutput

filedir=cd;
% names=getFileList(filedir,'tuning curves_OS data_All_cells',0,'anywhere');
names=getFileList(filedir,'tuning curves_Aug_27_2014',0,'anywhere');
if ~isempty(names)
    load(num2str(names{1}));
end

if ~isempty(VecX)
    matchind180=[];
    angleXYZ=rad2deg(real(acos(VecX.*VecXcomp+VecY.*VecYcomp+VecZ.*VecZcomp)));
    
% calculate mean output based on tuning curve
    for m=1:size(angleXYZ,2)
%         if ~isnan(angleXYZ(1,m))&& angleXYZ(1,m)<=90
        if ~isnan(angleXYZ(1,m))
            %output(m,1)=1-tuning.meanR(tuning.meanx==abs(floor(angleXYZ(1,m))-180)); % for old tuning curve
            output(m,1)=tuning.meanR(tuning.meanx==abs(floor(angleXYZ(1,m)))); % for new tuning curve
%         elseif ~isnan(angleXYZ(1,m))&& angleXYZ(1,m)>90
%              %output(m,1)=1-tuning.meanR(tuning.meanx==floor(angleXYZ(1,m))); % for old tuning curve
%             output(m,1)=tuning.meanR(tuning.meanx==floor(angleXYZ(1,m)-180)); % for new tuning curve
        else
            output(m,1)=nan;
        end
    end
    meanOutput=mean(output,'omitnan');
end    

%% calculate matchind matrices
if ~isempty(VecX)
    matchind180=[];
    matchindMix=[];
    matchindOSI=[];
    angleXYZ=rad2deg(real(acos(VecX.*VecXcomp+VecY.*VecYcomp+VecZ.*VecZcomp)));
    
    %vectorAngle(V1, V2);
    
    for c=10:10:10
        cutoff=(['c',num2str(c)]);
        matchind180.(cutoff)=100*(numel(find(angleXYZ<c))/numel(angleXYZ));
        matchindMix.(cutoff)=100*(((numel(intersect(find(angleXYZ>=90-c),find(angleXYZ<90+c)))/2)+numel(find(angleXYZ<c)))/numel(angleXYZ));
        %matchindMix.(cutoff)=100*(((sum(intersect(data.OSI(angleXYZ>=90-c),data.OSI(angleXYZ<90+c)))/2)+sum(data.OSI(angleXYZ<c)))/numel(angleXYZ)); 
        try matchindOSI.(cutoff)=100*(sum(data.OSI(angleXYZ<c))/numel(angleXYZ)); end
    end
    matchind90=[];
    angleXYZ90=[angleXYZ(find(angleXYZ<=90)),180-angleXYZ(find(angleXYZ>90))]; 
%     angleXYZ90 = vertcat(angleXYZ(find(angleXYZ <= 90)), 180 - angleXYZ(find(angleXYZ > 90)));
    
    %for c=5:5:45
    for c=10:10:10
        cutoff=(['c',num2str(c)]);
        matchind90.(cutoff)=100*(numel(find(angleXYZ90<c))/numel(angleXYZ90));
    end
end

% figure;
% subplot(1,2,1); histogram(angleXYZ,15);
% ax=gca; xlim([0 180]); ax.XTick=[0:45:180];
% subplot(1,2,2); histogram(angleXYZ90,15);
% ax=gca; xlim([0 90]); ax.XTick=[0:45:90];
% drawnow expose
%%  save 
allPossibleFields.field{k,h}=fieldID;
try allPossibleFields.anglediff180{k,h}=angleXYZ'; end
try allPossibleFields.anglediff90{k,h}=angleXYZ90'; end
allPossibleFields.VecXcomp{k,h}=VecXcomp';
allPossibleFields.VecYcomp{k,h}=VecYcomp';
allPossibleFields.VecZcomp{k,h}=VecZcomp';
try allPossibleFields.matchindOSI(k,h)=matchindOSI; end
try allPossibleFields.matchind180(k,h)=matchind180; end
try allPossibleFields.matchindMix(k,h)=matchindMix; end
try allPossibleFields.matchind90(k,h)=matchind90; end
try allPossibleFields.meanOutput(k,h)=meanOutput; end
if k==1 && h==1
    allPossibleFields.X=X';
    allPossibleFields.Y=Y';
    allPossibleFields.Z=Z';
    allPossibleFields.VecX=VecX';
    allPossibleFields.VecY=VecY';
    allPossibleFields.VecZ=VecZ';
end




end

