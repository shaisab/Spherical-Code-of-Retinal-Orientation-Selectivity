
%% calculates CI for the six coefficients and R2
realCoeff=xn180;
nullCoeff=null.coeff;

tails=2;                % tails=1 for one-tail and tails=2 for two-tails
alphavalue=0.01;        % alpha significance value
L=alphavalue/2;         % alpha value to be used for the lower confidence bound
U=1-(alphavalue/2);     % alpha value to be used for the upper confidence bound

% coefficients
qL=quantile(nullCoeff,L);       % qL and qU are the significance values for the L and U bounds, for the six coefficients
qU=quantile(nullCoeff,U); 
meanq=mean([qL;qU],1);
diff=meanq-qL;

% R square
qLR2=quantile(null.R2_180,L);   % qLR2 and qUR2 are the significance values for the L and U bounds, for the r square
qUR2=quantile(null.R2_180,U); 
meanqR2=mean([qLR2;qUR2],1);
diffR2=meanqR2-qLR2;

% plot coefficients and R2
figure;
ax = gca;
errorbar([meanq,meanqR2],[diff, diffR2],'k','LineStyle','none');
hold on
scatter(1:size(realCoeff,2)+1,[realCoeff,R2_180], 'r')
xlim([0.5 7.5]);
ylim([0 1]);
ax.XTickLabel = {'+ASC','+PSC','+LSC','-ASC','-PSC','-LSC','R^2'};


%% calculate rand-rand and real-rand matrix subtraction

matSize=size(null.randmat);
for r=2:matSize(2)
  randMinusRand(r-1)=max(max(null.randmat(1).mat180-null.randmat(r).mat180));
  realMinusRand(r-1)=max(max(mat180real-null.randmat(r).mat180));
end

figure;
plot(randMinusRand);
hold on
plot(realMinusRand);

%% calculate SD over many iteration of rand vectors compared to the real data
matSize=size(null.randmat);
for r=1:matSize(2)
  randSD(r)=std(std(null.randmat(r).mat180))/mean(mean(null.randmat(r).mat180));
  realSD(r)=std(std(mat180real))/mean(mean(mat180real));
end

figure;
histogram(randSD,15);
yl(1,1:2)=ylim;
hold on
plot(realSD,0:99,'r','LineWidth',3);
ylim([0 yl(:,2)]);
% xlabel('coefficient of variation');
% ylabel('number of itterations');

%% calculate the sum of all rand iterations
matSize=size(null.randmat);
randSum=zeros(size(null.randmat(1).mat180,1),size(null.randmat(1).mat180,2));
for r=1:matSize(2)
  randSum=randSum+null.randmat(r).mat180;
end
randAve=randSum./matSize(2);

alpha=0:resolution:360;
beta=0:resolution:180;
figure;
subplot(1,2,1)
contourf(alpha,beta,randAve);
subplot(1,2,2)
contourf(alpha,beta,mat180real);
subtitle('need to equalize the two colormaps');



