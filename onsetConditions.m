
speedterm=pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed;
a=(screenWidthum-widthum)/2;
b=roiLocumX;
c=radDendTreeum;
d=roiLocumY;
e=widthum;
f=heightum;
g=sqrt((f-d)^2+(b-e/2)^2);
beta=acosd((b-e/2)/g);
h=screenWidthum;
m=sqrt(h^2+h^2)/2;  % half diagonal of the screen




switch mergedmeta.metadataDB{1,1}.log.tested(p)
    case 0
        onset1=(a+b-c)/(speedterm);
        onset2=(a+b)/(speedterm);
    case 45
        x=g*sind(beta-45);
        n=m-x;
        onset1=(n-c)/(speedterm);
        onset2=n/(speedterm);
    case 180
        onset1=(a+(e-b)-c)/(speedterm);
        onset2=(a+(e-b))/(speedterm);
    case 225
        x=g*sind(beta-45);
        n=m+x;
        onset1=(n-c)/(speedterm);
        onset2=n/(speedterm);
    case 90
        onset1=(a+d-c)/(speedterm);
        onset2=(a+d)/(speedterm);
    case 135
        x=g*cosd(beta-45);
        n=m-x;
        onset1=(n-c)/(speedterm);
        onset2=n/(speedterm);
    case 270
        onset1=(a+(e-f)+(f-d)-c)/(speedterm);
        onset2=(a+(e-f)+(f-d))/(speedterm);
    case 315
        x=g*cosd(beta-45);
        n=m+x;
        onset1=(n-c)/(speedterm);
        onset2=n/(speedterm);
end



