
% Generate 50 random angles 
theta=randsample(360,50);       % the angles are given in degrees
rho=ones(50,1);                 % the radius equals to 1.

% create a polar histogram
f1=figure;
set(f1, 'color', [1 1 1]);
set(f1,'position',[50 50 500 500]);
polarhistogram(theta,40,'FaceColor',[0.29,0.66,0.03],'FaceAlpha',1,'EdgeColor',[0.30,0.75,0.93],'LineWidth',1);
%polarhistogram(theta,40,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',1,'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',1);
ax = gca;
ax.RTick = [];
ax.ThetaTick =[];
ax.ThetaColor = 'w';













