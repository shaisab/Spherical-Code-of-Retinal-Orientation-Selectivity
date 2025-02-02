% open the appropriate pooledmap files with Itrans

function [] = clusterBinnedCellsEmpiricalManual()

% change lobe to be either 'd' for dorsal lobe or 'n' for nasal lobe
lobe='d';


[Name,~] = uigetfile('*.mat','Select a pooledmap file with color-codded vectors');
            load(Name);

            
% % bins and bounds of bins for ON
% bins=[45,46,47,48,55,56,57,66,67,76]; 
% switch lobe
%     case 'd'
%         limbin=[80,140;60,120;60,110;30,90;65,130;60,130;60,120;70,130;60,120;40,90];   % dorsal lob
%     case 'n'
%         
% end


% bins and bounds of bins for ONOFF
bins=[27,28,45,46,47,48,55,56,57,66,67,76,77,78,87,36,37,38,33,43]; 
switch lobe
    case 'd'
        limbin=[30,90;30,90;80,140;60,120;60,110;30,90;65,130;60,130;60,120;70,130;60,120;40,90;60,90;60,120;60,90;50,100;40,100;30,90;120,150;100,140];  % dorsal lobe
    case 'n'
        limbin=[330,0;310,20;340,50;330,30;330,30;315,10;330,50;330,40;330,40;340,50;340,30;350,30;350,20;350,20;350,20;330,30;330,30;320,10;30,60;30,70];    % nasal lobe
end


for i=1:size(bins,2)
    bin=bins(i);
    [theta,~]=cart2pol(pooledmap.clusterdata{1,bin}.UVStvecComp.VecUSt,pooledmap.clusterdata{1,bin}.UVStvecComp.VecVSt);
    ang=mod(theta,2*pi);
    figure;
    circ_plot(ang);
    
    switch lobe
        case 'd'
            ind=logical(logical(rad2deg(ang)>limbin(i,1)).* logical(rad2deg(ang)<limbin(i,2)));
        case 'n'
            if limbin(i,1)>180
                ind=logical(logical(rad2deg(ang)>limbin(i,1))+ logical(rad2deg(ang)<limbin(i,2)));
            else
                ind=logical(logical(rad2deg(ang)>limbin(i,1)).* logical(rad2deg(ang)<limbin(i,2)));
            end
    end
    
    hold on;
    circ_plot(ang(ind),'r');
    meanAng(i)=mod(circ_mean(ang(ind)),2*pi);
    hold on;
    circ_plot(meanAng(i),'g'); 
    
    switch lobe
        case 'd'
            shift=meanAng(i)-pi/2;         % the shift needed to place the meanAng on the nasal lobe (0)
        case 'n'
            shift=meanAng(i)-0;      % the shift needed to place the meanAng on the dorsal lobe (pi/2)
    end
    
    angS{i}=ang-repmat(shift,size(ang,1),1);
    figure;
    circ_plot(angS{i});
    hold on;
    circ_plot(meanAng(i)-shift,'r'); 
  
end

%% plot a polar plot of angles
figure
hold on
for j=1:size(angS,2)
   circ_plot(angS{j},'k'); 
end

%% plot a quiver polar plot 

figure; 
hold on
for j=1:size(angS,2)
    [X,Y]=pol2cart(angS{j},1);
h=quiver(zeros(size(X,1),1),zeros(size(X,1),1),X,Y,1,'k');
h.ShowArrowHead='off';
end
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
axis equal;
    
end


