
d=deg2rad(0:45:315);
%tmpa=[10,30,40,80,34,23,20,10];
a=[10,80,10,10,10,80,10,10];
figure;plot(rad2deg(tmpd),tmpa)

[mu kappa] = circ_vmpar(alpha_rad);


 rad2deg(mu)
 
 
 R=Rmax*exp(k*cos((x-u)*pi/180))/exp(k);
 
 
 
 
 [THETA,RHO]=cart2pol(real(zm),imag(zm));
 
 rad2deg(THETA)
 
 
 % calculate pori according to Wang et al.,2010b
                        
 pori2=rad2deg(angle(sum(a.*exp(2i*deg2rad(d)))/sum(a)));
 
 
 
 ai=a(5:8)
 amean=mean([a(1:4);a(5:8)],1)
 di=d(1:4)
 pori=mean(di.*amean)/sum(amean)
 
 pori2=rad2deg(angle(sum(amean.*exp(2i*deg2rad(di)))/sum(amean)))
 