% pori_on_and_off

pori2=rad2deg(angle(sum(a.*exp(2i*deg2rad(d)))/sum(a)));
if pori2<=0
    pori2 = pori2+180;
end

oriFWHM,u,gof2]= doubleGaussianFit(ROIn,d,a);
pori=u;
oriR2=gof2.rsquare;