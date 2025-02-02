function [DX] = DiffX(m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function outputs the matrix we use for differentiating a surface
%   with respect to x. The input m is simply the size of the matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



B(:,1)=[ones(m-2,1)',0,0]';
B(:,2)=[-8*ones(m-1,1)',0]';
B(:,3)=[0,8*ones(m-1,1)']';
B(:,4)=[0,0,-1*ones(m-2,1)'];

d = [-2,-1,1,2];

DX = spdiags(B,d,m,m);

DX(1,1)=-25;
DX(1,2)=48;
DX(1,3)=-36;
DX(1,4)=16;
DX(1,5)=-3;

DX(2,1)=-3;
DX(2,2)=-10;
DX(2,3)=18;
DX(2,4)=-6;
DX(2,5)=1;

DX(m,m)=25;
DX(m,m-1)=-48;
DX(m,m-2)=36;
DX(m,m-3)=-16;
DX(m,m-4)=3;

DX(m-1,m)=3;
DX(m-1,m-1)=10;
DX(m-1,m-2)=-18;
DX(m-1,m-3)=6;
DX(m-1,m-4)=-1;

DX=1/(12)*DX;
 end

