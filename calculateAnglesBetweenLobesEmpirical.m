% open the appropriate pooledmap files with Itrans

function [] = calculateAnglesBetweenLobesEmpirical()

% input file and bins data
[Name,~] = uigetfile('*.mat','Select a pooledmap file with color-codded vectors');
            load(Name);
 
% % bins and bounds of bins for ON
% bins=[45,46,47,48,55,56,57,66,67,76]; 
% limbinD=[80,140;60,120;60,110;30,90;65,130;60,130;60,120;70,130;60,120;40,90];   % dorsal lob

% bins and bounds of bins for ONOFF
bins=[27,28,45,46,47,48,55,56,57,66,67,76,77,78,87,36,37,38,33,43]; 
limbinD=[30,90;30,90;80,140;60,120;60,110;30,90;65,130;60,130;60,120;70,130;60,120;40,90;60,90;60,120;60,90;50,100;40,100;30,90;120,150;100,140];  % dorsal lobe
limbinN=[330,0;310,20;340,50;330,30;330,30;315,10;330,50;330,40;330,40;340,50;340,30;350,30;350,20;350,20;350,20;330,30;330,30;320,10;30,60;30,70];    % nasal lobe

% find the cells included in each lobe
for i=1:size(bins,2)
    bin=bins(i);
    [theta,~]=cart2pol(pooledmap.clusterdata{1,bin}.UVStvecComp.VecUSt,pooledmap.clusterdata{1,bin}.UVStvecComp.VecVSt);
    ang=mod(theta,2*pi);
    figure;
    hold on;
    circ_plot(ang);
    indD=logical(logical(rad2deg(ang)>limbinD(i,1)).* logical(rad2deg(ang)<limbinD(i,2)));
    if limbinN(i,1)>180
        indN=logical(logical(rad2deg(ang)>limbinN(i,1))+ logical(rad2deg(ang)<limbinN(i,2)));
    else
        indN=logical(logical(rad2deg(ang)>limbinN(i,1)).* logical(rad2deg(ang)<limbinN(i,2)));
    end


% calculate mean angle of lobes
meanAngD(i)=mod(circ_mean(ang(indD)),2*pi);
meanAngN(i)=mod(circ_mean(ang(indN)),2*pi);

% calculate angle differences    
angDiff(i)=abs(meanAngN(i)-meanAngD(i));
if angDiff(i)>pi
    angDiff(i)=pi*2-angDiff(i);
end 

circ_plot(ang(indD),'r');
circ_plot(ang(indN),'g');
circ_plot(meanAngD(i),'k'); 
circ_plot(meanAngN(i),'k'); 
    
end    
    
[X,Y]=meshgrid(linspace(-1,1,10), linspace(-1,1,10));
X=reshape(X,100,1);
Y=reshape(Y,100,1);
x=X(bins);
y=Y(bins);

figure;
hold on;
angDiff=rad2deg(angDiff);
scatter(x,y,300,angDiff,'filled');
scatter(0,0,100,'k');
xlim([-1 1]); ylim([-1 1]);
axis equal
axis off
caxis([0 180]);

end


