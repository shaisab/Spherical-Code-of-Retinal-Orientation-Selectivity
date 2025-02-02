function [oriFWHM,u,k,gof2,x,xq,Rq]= doubleGaussianFitOSI(j,d,a)

x=d';
y=a';
% y_temp=a';
% y=y_temp./max(y_temp);

% for lower, upper, and startpoint the order of variables is: Rmax (maximum response), k (width of tuning curve), u (prederred orientation)
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Algorithm','Trust-Region',...
               'Display','iter',...
               'TolFun',1.0E-12,... 
               'TolX',1.0E-12,...
               'Lower',[max(a)-min(a),0,0],...
               'Upper',[3*max(a),135,180],...
               'StartPoint',[max(a)-min(a),45,90]);
%ft = fittype('(Rmax*(exp(k*cos((x-u)*pi/180))+exp(k*cos((x-u+180)*pi/180))))/exp(k)','options',fo);
ft = fittype('Rmax*(exp(k*cos((x-u)*pi/180))/exp(k))+Rmax*(exp(k*cos((x-u+180)*pi/180))/exp(k))','options',fo); 
% this function assumes: (1) the two peaks are 180 degrees apart, and (2)
% the peaks are of the same height. Both assumptions are reasonable for OS
% cells.

[curve2,gof2]=fit(x,y,ft)

% dispfigh2=figure(2000+j);clf;
% screendims=get(0,'MonitorPositions');
% if numel(screendims)==8
%     set(dispfigh2,'position',[1250 300 500 450]);
% elseif screendims(4)==800
%     set(dispfigh2,'position',[300 50 650 650]);
% else
%     set(dispfigh2,'position',[300 50 750 750]);
% end


% plot(x,y,'o');
% hold on
% plot(curve2,'m');

Rmax=curve2.Rmax;
k=curve2.k;
u=curve2.u;
RmaxC=Rmax*(1+exp(-2*curve2.k));

xq=0:315;
Rq = feval(curve2, xq);
% plot(xq, Rq,'k','LineWidth',2,'Color',[0 0 0]);


R=Rmax*(exp(k*cos((x-u)*pi/180))/exp(k))+Rmax*(exp(k*cos((x-u+180)*pi/180))/exp(k));

Rq = interp1(x,R,xq,'spline');
Rq=Rq-min(R);
R=mat2gray(R);
Rq=mat2gray(Rq);

% plot(x,R,'o',xq,Rq,':.');
% Rq_linear=interp1(x,R,xq,'linear');
% plot(x,R,'o',xq,Rq_linear);
% hold off
% drawnow expose

ind=find(Rq(90:225)==max(Rq(90:225)));      % find a maximum point around the middle of the curve
RqTmp=Rq(ind+89-89:ind+89+89);              % truncate the curve such that it shows only one peak
xqTmp=xq(ind+89-89:ind+89+89);

% RqTmp=Rq(ind+89-89:ind+89+89);              % truncate the curve such that it shows only one peak
% xqTmp=xq(ind+89-89:ind+89+89);

% ind=find(Rq(90:300)==max(Rq(90:300)));      % find a maximum point around the middle of the curve
% RqTmp=Rq(ind-50:ind+100);              % truncate the curve such that it shows only one peak
% xqTmp=xq(ind-50:ind+100);
oriFWHM=fwhmOSI(xqTmp,RqTmp');                  % calculate FWHM of the function in degrees

%% code to visualize the function
% Rmax=1.7;
% k=90;
% u=180;
% b=min(a);
% R=(Rmax*(exp(k*cos((x-u)*pi/180))+exp(k*cos((x-u+180)*pi/180))))/exp(k);
% %R=Rmax*exp(k*cos((x-u)*pi/180))/exp(k);
% figure;
% plot(x,R,'-o');


end