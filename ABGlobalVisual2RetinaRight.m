function [alpha,beta,Normal]= ABGlobalVisual2RetinaRight(Normal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This code converts a normal passing from the ORIGIN OUTWARD in global 
%   visual spaceto a normal passing from the ORIGIN OUTWARD in the retinal 
%   Space and then converts this to a normal and alpha beta on the retina. 
%
%   Inputs
%   1. Normal - This is a normal measured from the origin outward in the
%   global visual space. NOTE: Normal must be passed as a column vector.
%
%   Outputs
%   1. alpha - The alpha coordinate as measured on the retina.
%   2. beta - the beta cordinate as measured on the retina.
%   3. normal - This is the normal relative to the retina.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rotation matrix 

% The rotation matrix is computed in Mathematica notebook: RetinaNotebook.
% The code in Mathematica computes the rotation from a retina with optical
% axis [0,0,1] into the optical axis in global visual space and then
% inverts this matrix. The sequence of rotations is the following:
%
%   1. Rotate about the opitical axis ([0,0,1]) by 90 degrees to align the
%   muscle with the x-axis (lambda-bregma). This is necessary since we 
%   assume the cuts are aligned with the x-axis in the Cartesian coordinate 
%   system about the retina.
%
%   2. Rotate by 35.6 degrees about the opitcal axis by 35.6 degrees to
%   account for the fact that the muscles are 35.6 degrees from the x-axis.
%
%   3. Rotate by -68 degrees about the y axis. This rotation will align the
%   optical axis to be 22 degrees above the horizon. 
%
%   4. Rotate by -64 degrees about the z-axis. This rotation will align the
%   optcal axis to be 64 degrees away from the mid sagital plane. 

Rot=[0.907945,0.414229,-0.0636356;...
-0.102142,0.365987,0.924998;...
0.406451,-0.833347,0.374607];           % rot while acounting for 2.5 deg correction

Normal=Rot*Normal;

%% Conversion to alpha beta coordinates. 
% This calculation is done based off of the following parametrization of
% the sphere:
%
%   x=cos(alpha)sin(beta)
%   y=sin(alpha)sin(beta)
%   z=-cos(beta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=atan2(Normal(2),Normal(1));
beta=acos(-Normal(3));

% correction for rotation of the eye in the orbit based on the angle
% difference Lambda-Bregma axis and the medial and lateral muscles measured
% by me and Ananya, file excel name: 'extra ocular muscles vs
% Lambda-Bregma'
% in directory: C:\Users\ssabbah\Documents\NathanAnatomy\extra ocular muscles vs Lambda-Bregma
%correction=deg2rad(-35.62);
%alpha=alpha+correction;

%alpha=alpha+deg2rad(2.5);        % adding the median alpha correction


alpha=rad2deg(mod(alpha,2*pi));
beta=rad2deg(mod(beta,pi));


% standardize alpha to range between 0 and 360 degrees, and beta to range between 90-180 degrees
% if rad2deg(beta)>90
%     beta=rad2deg(beta);
%     if rad2deg(alpha)>=0
%         alpha=rad2deg(alpha);
%     else
%         alpha=360+rad2deg(alpha);
%     end    
% else
%     beta=180-rad2deg(beta);
%     if rad2deg(alpha)<=180
%         alpha=180+rad2deg(alpha);
%     else rad2deg(alpha)>180 
%         alpha=rad2deg(alpha)-180;
%     end
% end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%