
%% for ON-OFF DSGCs
figure
y=[0.92,0.97,0.97;...
0.42,0.75,0.8;...
0.37,0.74,0.82];
%y = [2 4 6; 3 4 5];
b = bar(y);
%labels={'Translation','Rotation','Hybrid'};
labels={'Whole','Temporal','Extreme temporal'};
set(gca, 'XTick', 1:3, 'XTickLabel', labels,'FontSize',14);
xlabel('Retinal region','FontSize',20);
ylabel('\itR\rm^2','FontSize',20);

%% for ON DSGCs
figure
y=[0.87,0.95,0.96;...
0.29,0.79,0.86;...
0.33,0.53,0.58];
%y = [2 4 6; 3 4 5];
b = bar(y);
%labels={'Translation','Rotation','Hybrid'};
labels={'Whole','Temporal','Extreme temporal'};
set(gca, 'XTick', 1:3, 'XTickLabel', labels,'FontSize',14);
xlabel('Retinal region','FontSize',20);
ylabel('\itR\rm^2','FontSize',20);

