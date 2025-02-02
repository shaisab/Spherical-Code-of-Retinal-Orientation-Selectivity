

function F = myfunSurf6c(v,C1,C2,C3,C4,C5,C6,Cr)
a1=v(1);
a2=v(2);
a3=v(3);
a4=v(4);
a5=v(5);
a6=v(6);

F=abs(Cr-(a1.*C1+a2.*C2+a3.*C3+a4.*C4+a5.*C5+a6.*C6));

end
