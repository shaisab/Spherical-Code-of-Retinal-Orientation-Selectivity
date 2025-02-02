function [E] = SectorEnergy(input,DS,DTH,nu,S,s,R,n,ds,dth,a,b,a1,a2)

%% Conversion of input into functions to be varied over
rho=input(:,1:n);
f=input(:,n+1:2*n);

rho(1,1:a)=s(1:a);
rho(n,1:b)=s(1:b);

f(1,1:a)=a1;
f(n,1:b)=a2;

DTrho=DTH*rho;
DSrho=(DS*rho')';

DTf=DTH*f;
DSf=(DS*f')';

%% Definition of strain terms.
gamma11=DSrho.^2+rho.^2.*DSf.^2-1;
gamma12=(DSrho.*DTrho+rho.^2.*DTf.*DSf)./(R*sin(S/R));
gamma22=(DTrho.^2+rho.^2.*DTf.^2-R^2.*sin(S/R).^2)./(R*sin(S/R)).^2;

E(:,1:n)=sqrt(nu)*(gamma11+gamma22)*sqrt(ds*dth).*(R*sin(S/R)).^(1/2);
E(:,n+1:2*n)=sqrt((1-nu))*gamma11*sqrt(ds*dth).*(R*sin(S/R)).^(1/2);
E(:,2*n+1:3*n)=sqrt((1-nu))*gamma22*sqrt(ds*dth).*(R*sin(S/R)).^(1/2);
E(:,3*n+1:4*n)=sqrt(2)*sqrt(1-nu)*gamma12*sqrt(ds*dth).*(R*sin(S/R)).^(1/2);


E(1,1:4*n)=E(1,1:4*n)/sqrt(2);
E(n,1:4*n)=E(n,1:4*n)/sqrt(2);

E(:,1)=E(:,1)/sqrt(2);
E(:,n)=E(:,n)/sqrt(2);

E(:,n+1)=E(:,n+1)/sqrt(2);
E(:,2*n)=E(:,2*n)/sqrt(2);

E(:,2*n+1)=E(:,2*n+1)/sqrt(2);
E(:,3*n)=E(:,3*n)/sqrt(2);

E(:,3*n+1)=E(:,3*n+1)/sqrt(2);
E(:,4*n)=E(:,4*n)/sqrt(2);



end

