


pT1=[-0/sqrt(1.^2+0.^2),1/sqrt(1.^2+0.^2),0];     % a point on the plane perpendicular to translation axis 1

pOpt=[-1/sqrt(0.^2+1.^2),0/sqrt(0.^2+1.^2),0];     % a point on the plane perpendicular to the optic axis

[P,N,check]=plane_intersect(RightoptL,pOpt,[channels.xT1(end),channels.yT1(end),channels.zT1(end)],pT1)
