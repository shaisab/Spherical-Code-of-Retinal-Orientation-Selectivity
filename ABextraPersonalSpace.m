function [alpha,beta]= ABextraPersonalSpace(Normal)

%OptAxis=[0.406,-0.833,0.375];       % optL
%OptAxis=[0.415,-0.7215,0.553];     % optNr


% Rot=vrrotvec2mat(vrrotvec(OptAxis',[0,0,1]'));
% 
% Normal=Rot*Normal';

alpha=atan2(Normal(2),Normal(1));
beta=acos(-Normal(3));

% correction for rotation of the eye in the orbit based on the angle
% difference Lambda-Bregma axis and the medial and lateral muscles measured
% by me and Ananya, file excel name: 'extra ocular muscles vs
% Lambda-Bregma'
% in directory: C:\Users\ssabbah\Documents\NathanAnatomy\extra ocular muscles vs Lambda-Bregma
% correction=deg2rad(-35.62);
% alpha=alpha+correction;

alpha=rad2deg(alpha)
beta=rad2deg(beta)

% standardize alpha to range between 0 and 360 degrees, and beta to range between 90-180 degrees
if beta>90
    beta=180-beta;
    if alpha>=0
        alpha=360-alpha;
    else
        alpha=180+alpha;
    end    
else
    beta=180-beta;
    if alpha<=180
        alpha=-alpha;
    else alpha>180 
        alpha=alpha-180;
    end
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%