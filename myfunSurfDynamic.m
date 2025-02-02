

function F = myfunSurfDynamic(v,regressor,Cr)

a=v;
% a1=v(1);
% a2=v(2);
% a3=v(3);
% a4=v(4);
% a5=v(5);
% a6=v(6);

term=0;
for k=1:size(v,2)
term=term+a(k).*regressor(:,k);
end

%F=abs(Cr-(a1.*C1+a2.*C2+a3.*C3+a4.*C4+a5.*C5+a6.*C6));
F=abs(Cr-(term));
end
