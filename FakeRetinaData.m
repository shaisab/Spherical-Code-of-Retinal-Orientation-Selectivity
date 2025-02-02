
%% create fake data (position only) uniformly across the retina
phys=dlmread('phys');
R=phys(1);
M=phys(2); %This should be normalized to 1.

s=linspace(0,M,45);
th=linspace(0,360,45);

th=deg2rad(th);

k=0;
for i=1:45,
    for j=1:45,
        k=k+1;
   x(k)=R*cos(th(j))*sin(s(i)/R);
   y(k)=R*sin(th(j))*sin(s(i)/R);
   z(k)=R-R*cos(s(i)/R);
    end
   
end

plot3(x,y,z,'.');
axis equal

pooledmap.XYZvecComp.X=x';
pooledmap.XYZvecComp.Y=y';
pooledmap.XYZvecComp.Z=z';
pooledmap.XYZvecComp.VecX=[];
pooledmap.XYZvecComp.VecY=[];
pooledmap.XYZvecComp.VecZ=[];

%% save pooledmap file
filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
discPoints=40;
save(['pooledMap_',num2str(filename),'_Fake data_',num2str(discPoints),'_created_',currenttime,'.mat'],'pooledmap');
    

