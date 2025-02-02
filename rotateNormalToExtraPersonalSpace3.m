function [alphaRG,betaRG]=rotateNormalToExtraPersonalSpace3(side,alpha,beta)

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

%% find the rotation matrix needed to transform the right (RightoptL) or left (LeftoptL) optic axis  to the modified LB, and use it to rotate all vectors
% calculate the alpha beta on the left retina based on the alpha beta in
% the right retina
localOpt=[0,0,1];
Hor=[1,0,0];

hold on;
quiver3(0,0,0,localOpt(1),localOpt(2),localOpt(3),'r','LineWidth',4);
quiver3(0,0,0,Hor(1),Hor(2),Hor(3),'k','LineWidth',6);

%acosd(dot(Normal,localOpt))
% optic axes are defined relative to the horizon in the ambulatory mouse with lambda-bregma tilted 29 deg down.
RightoptL=[0.406450649610760,-0.833347328309341,0.374606593415912];     
LeftoptL=[0.406450649610760,0.833347328309341,0.374606593415912];

switch side
    case 'Right'
        % convert the normal in a right eye system into a normal in  the
        % global system
        R3=vrrotvec2mat(vrrotvec([0,0,1]',RightoptL'));
        Normalr=R3*Normal';
    case 'Left'
        % convert the normal in a left eye system into a normal in  the
        % global system
        R3=vrrotvec2mat(vrrotvec([0,0,1]',LeftoptL'));
        Normalr=R3*Normal';
        %Normalr=-Normalr;
end
localOptr=R3*localOpt';
Horr=R3*Hor';

quiver3(0,0,0,localOptr(1),localOptr(2),localOptr(3),'b','LineWidth',4);
quiver3(0,0,0,Horr(1),Horr(2),Horr(3),'k','LineWidth',4);
        
%% correction to rotate the Horr to be horizontal again. The rotation to the optic was not correct, resulting in hor not being horizontal anymore.
% to correct for this error, we measured the angle between the new Hor
% after rotation (Horr) and the true horizontal (1,0,0). We calculated the
% rotation matrix needed to rotate the Horr to be horizontal again. We
% applied the same rotation matrix to rotate the normal.

x=linspace(-1,1,30);
y=linspace(-1,1,30);
[X,Y]=meshgrid(x,y);
Z=zeros(30,30);

hold on
surf(X,Y,Z);
surf(X,Y,-1.08516*X+2.22371*Y);
axis square;
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
xlabel('x');ylabel('y');zlabel('z');

switch side
    case 'Right'
        Ref=[0.8987;0.4385;0];
    case 'Left'
        Ref=[-0.8987;0.4385;0];
end
Horr_Ref=acosd(dot(Horr,Ref)/(norm(Horr)*norm(Ref)));
quiver3(0,0,0,Ref(1),Ref(2),Ref(3),'g','LineWidth',4)

% HorrP=[Horr(1),Horr(2),0];
% Horr_Ref=acosd(dot(Horr,HorrP)/(norm(Horr)*norm(HorrP)));
% quiver3(0,0,0,HorrP(1),HorrP(2),HorrP(3),'g','LineWidth',4)

rot=AxelRot(Horr_Ref,localOptr,[0,0,0]);
rot(:,4)=[];
rot(4,:)=[];

Horr2=rot*Horr;
localOptr=rot*localOptr;

quiver3(0,0,0,Normalr(1),Normalr(2),Normalr(3),'y','LineWidth',4)
Normalr=rot*Normalr;

quiver3(0,0,0,Horr2(1),Horr2(2),Horr2(3),'--k','LineWidth',4)
quiver3(0,0,0,Normalr(1),Normalr(2),Normalr(3),'--y','LineWidth',4)
quiver3(0,0,0,localOptr(1),localOptr(2),localOptr(3),'--m','LineWidth',4)

%% correction for rotation of the eye in the orbit based on the angle
% difference Lambda-Bregma axis and the medial and lateral muscles measured
% by me and Ananya, file excel name: 'extra ocular muscles vs
% Lambda-Bregma'
% in directory: C:\Users\ssabbah\Documents\NathanAnatomy\extra ocular muscles vs Lambda-Bregma
correction=35.62;     %;

switch side
    case 'Right'
        rot=AxelRot(correction,localOptr,[0,0,0]);
        rot(:,4)=[];
        rot(4,:)=[];
    case 'Left'
        rot=AxelRot(-correction,localOptr,[0,0,0]);
        rot(:,4)=[];
        rot(4,:)=[];  
end

Normalr=rot*Normalr;
switch side
    case 'Right'
        [alpha,beta]=ABextraPersonalSpaceRight(Normalr);
    case 'Left'
        [alpha,beta]=ABextraPersonalSpaceLeft(Normalr);
end
alphaRG=alpha;
betaRG=beta;

end