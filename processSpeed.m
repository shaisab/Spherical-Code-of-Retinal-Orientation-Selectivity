

numDS=[sum(speed015(:,1)),sum(speed03(:,1)),sum(speed06(:,1)),sum(speed09(:,1))];
numOS=[sum(speed015(:,2)),sum(speed03(:,2)),sum(speed06(:,2)),sum(speed09(:,2))];

% calculate the DSI of all DS cells at speed 0.6
DSI(:,1)=speed015((speed06(:,1)==1),4);
DSI(:,2)=speed03((speed06(:,1)==1),4);
DSI(:,3)=speed06((speed06(:,1)==1),4);
DSI(:,4)=speed09((speed06(:,1)==1),4);
[~,maxDSI]=max(DSI,[],2);
maxDSI015=sum(maxDSI==1);
maxDSI03=sum(maxDSI==2);
maxDSI06=sum(maxDSI==3);
maxDSI09=sum(maxDSI==4);

j=1;
for i=1:size(speed015,1)
 if (speed015(i,1)==1) || (speed03(i,1)==1) || (speed06(i,1)==1) || (speed09(i,1)==1)
    sDS(j,:)=[speed015(i,1),speed03(i,1),speed06(i,1),speed09(i,1),i];
    j=j+1;
 end
end

speed09((speed06(:,1)<speed09(:,1)),3)      % the ONorONOFForOFF state of cells that are DS with speed09 but not with speed06
speed03((speed06(:,1)<speed03(:,1)),3)      % the ONorONOFForOFF state of cells that are DS with speed03 but not with speed06
speed015((speed06(:,1)<speed015(:,1)),3)    % the ONorONOFForOFF state of cells that are DS with speed015 but not with speed06


j=1;
for i=1:size(speed015,1)
 if (speed015(i,2)==1) || (speed03(i,2)==1) || (speed06(i,2)==1) || (speed09(i,2)==1)
    sOS(j,:)=[speed015(i,2),speed03(i,2),speed06(i,2),speed09(i,2),i];
    j=j+1;
 end
end



speed06(speed06(:,2).*speed06(:,3)),3)  

