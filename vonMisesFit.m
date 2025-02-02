function []= vonMisesFit(d,a)

%d=(0:45:315);
%a=[10,80,10,10,10,80,10,10];
%amean=mean([a(1:4);a(5:8)],1);
%amean=[amean,amean(1)];
%dmean=d(1:5);
%dmean=d(1:4);

%Rmax=max(amean);

% for lower, upper, and startpoint the order of variables is: Rmax (maximum response), b
% (baseline), k (width of tuning curve), u (prederred orientation)
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Algorithm','Trust-Region',...
               'Display','iter',...
               'TolFun',1.0E-20,... 
               'TolX',1.0E-20,...
               'Lower',[0,0,0,0],...
               'Upper',[3*max(a),min(a),90,180],...
               'StartPoint',[0,min(a),180,90]);
%ft = fittype('(Rmax*exp(k*cos((x-u)*pi/180))/exp(k))','problem','Rmax','options',fo);
%ft = fittype('(Rmax*exp(k*cos((x-u)*pi/180))/exp(k))','options',fo);
%ft = fittype('(b+Rmax*(exp(k*cos((x-u)*pi/180))+exp(k*cos((x-u+180)*pi/180))))/exp(k)','problem','b','options',fo);
%ft = fittype('(Rmax*exp(k*cos((x-u)*pi/180))/exp(k))+(Rmax*exp(k*cos((x-u+180)*pi/180))/exp(k))','options',fo);
ft = fittype('b+(Rmax*(exp(k*cos((x-u)*pi/180))+exp(k*cos((x-u+180)*pi/180))))/exp(k)','options',fo);

% x=dmean';
% y=amean';

x=d';
y=a';

%[curve2,gof2]=fit(x,y,ft,'problem',max(y))
%[curve2,gof2]=fit(x,y,ft,'problem',min(y))
[curve2,gof2]=fit(x,y,ft)

figure;
plot(x,y,'o');
hold on
plot(curve2,'m');
hold off

%% code to visualize the function
% Rmax=1.7;
% k=90;
% u=180;
% b=min(a);
% R=b+Rmax*(exp(k*cos((x-u)*pi/180))+exp(k*cos((x-u+180)*pi/180)))/exp(k);
% %R=Rmax*exp(k*cos((x-u)*pi/180))/exp(k);
% figure;
% plot(x,R,'-o');


end