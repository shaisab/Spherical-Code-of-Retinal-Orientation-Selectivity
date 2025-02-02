function [alphaG,betaG,Normalr]=rotateNormalToExtraPersonalSpace4(side,alpha,beta)

% side can be either 'Right' or 'Left'

% invert the normal to run from the singularity on the retina toward the
% nodal point
% alpha=mod(alpha+180,360);
% beta=180-beta;

%% tarnsform alpha of the left eye to right-hand system (because for example, clockwise rotation on the right eye would be counterclockwise rotation on the left eye)
if strcmp(side,'Left')
alpha=-alpha;
alpha=mod(alpha,360);
end
%% calculate the normal based on alpha and beta, while tarnsforming beta of the left eye to right-hand system
N1=cosd(alpha)*sind(beta);
N2=sind(alpha)*sind(beta);
N3=-cosd(beta);

Normal=[N1,N2,N3];

normV=norm(Normal);
Normal(1)=Normal(1)/normV;
Normal(2)=Normal(2)/normV;
Normal(3)=Normal(3)/normV;

%% rotation Right
% rotR and RotL were computed in Mathematica retina notebook
switch side
    case 'Right'
         rotR=[0.907945,-0.102142,0.406451;...
        0.414229,0.365987,-0.833347;...
        -0.0636356,0.924998,0.374607];
        Normalr=rotR*Normal';
    case 'Left'
        rotL=[0.907945,0.102142,0.406451;...
        -0.414229,0.365987,0.833347;...
        -0.0636356,-0.924998,0.374607];
        Normalr=rotL*Normal';
end

switch side
    case 'Right'
        [alphaG,betaG]=ABextraPersonalSpaceRight(Normalr);
    case 'Left'
        [alphaG,betaG]=ABextraPersonalSpaceLeft(Normalr);
end

end