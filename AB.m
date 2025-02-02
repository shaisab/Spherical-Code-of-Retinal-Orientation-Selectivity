function [alpha,beta]= AB(Normal)

      
OptAxis=[0.406450649610760,-0.833347328309341,0.374606593415912];    % optL from Oct_28

Rot=vrrotvec2mat(vrrotvec(OptAxis',[0,0,1]'));

Normal=Rot*Normal';

alpha=atan2(Normal(2),Normal(1));
beta=acos(-Normal(3));

% correction for rotation of the eye in the orbit based on the angle
% difference Lambda-Bregma axis and the medial and lateral muscles measured
% by me and Ananya, file excel name: 'extra ocular muscles vs
% Lambda-Bregma'
% in directory: C:\Users\ssabbah\Documents\NathanAnatomy\extra ocular muscles vs Lambda-Bregma
correction=deg2rad(-35.62);
alpha=alpha+correction;

alpha=alpha+deg2rad(2.5);        % adding the median alpha correction


% standardize alpha to range between 0 and 360 degrees, and beta to range between 90-180 degrees
if rad2deg(beta)>90
    beta=rad2deg(beta);
    if rad2deg(alpha)>=0
        alpha=rad2deg(alpha);
    else
        alpha=360+rad2deg(alpha);
    end    
else
    beta=180-rad2deg(beta);
    if rad2deg(alpha)<=180
        alpha=180+rad2deg(alpha);
    else rad2deg(alpha)>180 
        alpha=rad2deg(alpha)-180;
    end
%     if rad2deg(alpha)<=180 && rad2deg(alpha)>=0 
%         alpha=180+rad2deg(alpha);
%     elseif rad2deg(alpha)>180 
%         alpha=rad2deg(alpha)-180;
%     elseif rad2deg(alpha)<0
%         alpha=360+rad2deg(alpha); 
%     end
end

alpha=deg2rad(alpha);
beta=deg2rad(beta);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%