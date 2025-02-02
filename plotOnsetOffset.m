function [] = plotOnsetOffset(physdata,k)

if strcmp(physdata.log(1).expType,'DB')
    screenWidthum=3000;  % 3 mm at the plane of the retina
    pixelprojwidth=screenWidthum/600;   % 3 mm / 600 pixels
    radDendTreeum=300;
    widthum=256;         % width of FOV in microns
    roiLocumX=128;       % roiLocumX and roiLocumY assume cell is located in the middle of the FOV 
    roiLocumY=64;
    corrSpeedFactor=1.1;
    
    switch (k-1)*45         % loop on all 8 directions
        case 0
            onset1=(((screenWidthum-widthum)/2)+roiLocumX-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+roiLocumX)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
        case 45
            onset1=(((screenWidthum-widthum)/2)+roiLocumX-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+roiLocumX)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
        case 180
            onset1=(((screenWidthum-widthum)/2)+widthum-roiLocumX-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+widthum-roiLocumX)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
        case 225
            onset1=(((screenWidthum-widthum)/2)+widthum-roiLocumX-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+widthum-roiLocumX)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
        case 90
            onset1=(((screenWidthum-widthum)/2)+roiLocumY-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+roiLocumY)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
        case 135
            onset1=(((screenWidthum-widthum)/2)+roiLocumY-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+roiLocumY)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
        case 270
            onset1=(((screenWidthum-widthum)/2)+widthum-roiLocumY-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+widthum-roiLocumY)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
        case 315
            onset1=(((screenWidthum-widthum)/2)+widthum-roiLocumY-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
            onset2=(((screenWidthum-widthum)/2)+widthum-roiLocumY)/(pixelprojwidth*corrSpeedFactor*60*physdata.log(1).speed);
    end
    offset1=onset1+(physdata.log(1).width/(corrSpeedFactor*60*physdata.log(1).speed));
    offset2=onset2+(physdata.log(1).width/(corrSpeedFactor*60*physdata.log(1).speed));
    
    onset1=onset1+physdata.timeAdd(1);
    onset2=onset2+physdata.timeAdd(1);
    offset1=offset1+physdata.timeAdd(1);
    offset2=offset2+physdata.timeAdd(1);
    
    
    
    hold on;
    plot([onset1 onset1],[0 100],'m');
    hold on;
    plot([onset2 onset2],[0 100],'m');
    hold on;
    plot([offset1 offset1],[0 100],'m');
    hold on;
    plot([offset2 offset2],[0 100],'m');
end

end