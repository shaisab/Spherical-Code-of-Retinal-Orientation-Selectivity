function [alpha,beta]= ABextraPersonalSpaceRight(Normal)

alpha=atan2(Normal(2),Normal(1));
beta=acos(-Normal(3));

beta=rad2deg(beta);
alpha=rad2deg(alpha);
alpha=mod(alpha,360);



end
