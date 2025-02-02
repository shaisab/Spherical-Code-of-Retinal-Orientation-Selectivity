

function F = myfunSurf(v,C1,C2,C3,Cr)
a1=v(1);
a2=v(2);
a3=v(3);

F=abs(Cr-(a1.*C1+a2.*C2+a3.*C3));

end
