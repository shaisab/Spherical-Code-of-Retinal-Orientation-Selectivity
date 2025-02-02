function [dirFWHM,u,k,gof2,x,R,xq,Rq]= singleGaussianFit(ROIn,d,a)

x=d';
y=a';

% for lower, upper, and startpoint the order of variables is: Rmax (maximum response), k (width of tuning curve), u (prederred orientation)
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Algorithm','Trust-Region',...
               'Display','iter',...
               'TolFun',1.0E-12,... 
               'TolX',1.0E-12,...
               'Lower',[max(a)-min(a),0,0],...
               'Upper',[3*max(a),360,180],...
               'StartPoint',[max(a)-min(a),45,90]);
ft = fittype('Rmax*(exp(k*cos((x-u)*pi/180))/exp(k))','options',fo); 
% this function assumes a single peak, reasonable for DS cells
% cells.

[curve2,gof2]=fit(x,y,ft)
%plot(x,y,'o');
hold on
%plot(curve2,'m');

Rmax=curve2.Rmax;
k=curve2.k;
u=curve2.u;
RmaxC=Rmax*(1+exp(-2*curve2.k));


R=Rmax*(exp(k*cos((x-u)*pi/180))/exp(k));
if size(x,1)<=8
xq=0:360;
x=[x;360];
R=[R;R(1)];
Rq=interp1(x,R,xq,'spline');
R=mat2gray(R);
Rq=mat2gray(Rq);
plot(x,R,'o',xq,Rq,'-');
drawnow expose
else
 R=mat2gray(R);
 xq=x;
 Rq=R;
 plot(x,R,'-');
end

dirFWHM=fwhm(xq,Rq);

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