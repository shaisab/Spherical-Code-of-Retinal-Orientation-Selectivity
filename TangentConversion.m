function [VOut] = TangentConversion(x,y,z,R,V)

%Projection 1
P1x=x/sqrt(x^2+y^2)*(R-z)/R;
P1y=x/sqrt(x^2+y^2)*(R-z)/R;
P1z=sqrt((2*R-z)*z)/R;

P1=[P1x,P1y,P1z];
P1=P1/norm(P1);

comp1=V.*P1;

%Projection 2
P2x=-y/sqrt(x^2+y^2)*sqrt((2*R-z)*z);
P2y=x/sqrt(x^2+y^2)*sqrt((2*R-z)*z);
P2z=0;

P2=[P2x,P2y,P2z];
P2=P2/norm(P2);

comp2=V.*P2;

VOut=comp1.*P1+comp2.*P2;



end

