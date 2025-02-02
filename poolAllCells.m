
function startAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to startAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global fieldType
global visAngle
global ONorONOFForOFFstate
global CARTstate
global RBPMSstate
global INJECTstate
global RETROstate
global BGstate
global COLORstate
global CLUSTERstate
global ALPHACORRstate
global OSIlow
global OSIhigh

discPoints=str2num(get(handles.discPoints,'String'));

% determine the current ONorONOFForOFFstate
ONorONOFForOFFstate=[];
if get(handles.ONOSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,1]; end      % ONDS
if get(handles.ONOFFOSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,2]; end   %ONOFFDS
if get(handles.OFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,3]; end     % OFFDS
if get(handles.ONtDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,4]; end     % ONtDS
if get(handles.ONsusOFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,5]; end     % ONsusOFFDS

% determine the current INJECTstate
INJECTstate=[];
if get(handles.MTNDbutton,'value')==1 INJECTstate=[INJECTstate,1]; end      % MTNd
if get(handles.SFbutton,'value')==1 INJECTstate=[INJECTstate,2]; end        % SF
if get(handles.NOTbutton,'value')==1 INJECTstate=[INJECTstate,3]; end       % NOT
if get(handles.SCbutton,'value')==1 INJECTstate=[INJECTstate,4]; end        % SC

%initializing variables for (X,Y,Z) coordinates
pooledmap.XYZvecComp.X=[];
pooledmap.XYZvecComp.Y=[];
pooledmap.XYZvecComp.Z=[];
pooledmap.XYZvecComp.VecX=[];
pooledmap.XYZvecComp.VecY=[];
pooledmap.XYZvecComp.VecZ=[];
pooledmap.XYZvecComp.VecXcompASC=[];
pooledmap.XYZvecComp.VecYcompASC=[];
pooledmap.XYZvecComp.VecZcompASC=[];
pooledmap.XYZvecComp.angleXYZASC=[];
pooledmap.XYZvecComp.VecXcompPSC=[];
pooledmap.XYZvecComp.VecYcompPSC=[];
pooledmap.XYZvecComp.VecZcompPSC=[];
pooledmap.XYZvecComp.angleXYZPSC=[];
pooledmap.XYZvecComp.VecXcompLSC=[];
pooledmap.XYZvecComp.VecYcompLSC=[];
pooledmap.XYZvecComp.VecZcompLSC=[];
pooledmap.XYZvecComp.angleXYZLSC=[];
%initializing variables for UV coordinates
pooledmap.UVvecComp.U=[];
pooledmap.UVvecComp.V=[];
pooledmap.UVvecComp.VecU=[];
pooledmap.UVvecComp.VecV=[];
pooledmap.UVvecComp.VecUcompASC=[];
pooledmap.UVvecComp.VecVcompASC=[];
pooledmap.UVvecComp.angleUVASC=[];
pooledmap.UVvecComp.VecUcompPSC=[];
pooledmap.UVvecComp.VecVcompPSC=[];
pooledmap.UVvecComp.angleUVPSC=[];
pooledmap.UVvecComp.VecUcompLSC=[];
pooledmap.UVvecComp.VecVcompLSC=[];
pooledmap.UVvecComp.angleUVLSC=[];
%initializing variables for UVstandard coordinates
pooledmap.UVStvecComp.USt=[];
pooledmap.UVStvecComp.VSt=[];
pooledmap.UVStvecComp.VecUSt=[];
pooledmap.UVStvecComp.VecVSt=[];
pooledmap.UVStvecComp.VecUStcompASC=[];
pooledmap.UVStvecComp.VecVStcompASC=[];
pooledmap.UVStvecComp.angleUVStASC=[];
pooledmap.UVStvecComp.VecUStcompPSC=[];
pooledmap.UVStvecComp.VecVStcompPSC=[];
pooledmap.UVStvecComp.angleUVStPSC=[];
pooledmap.UVStvecComp.VecUStcompLSC=[];
pooledmap.UVStvecComp.VecVStcompLSC=[];
pooledmap.UVStvecComp.angleUVStLSC=[];


i=0;        % the index in the pooledmap variables

filedir=cd;

switch ALPHACORRstate
    case 1
        namesASC=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_ASC'],0,'anywhere');
        
        for k=1:size(namesASC,2)
            ind=strfind(namesASC{k},'_');
            expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
            
            namesMap=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f')],0,'anywhere');
            for h=1:size(namesMap,2)
                if sum(strfind(namesMap{h},['ASC_',expDate{k}])) && sum(strfind(namesMap{h},'ASC'))
                    vecCompASC{k}=load(num2str(namesMap{h}));
                elseif sum(strfind(namesMap{h},['PSC_',expDate{k}])) && sum(strfind(namesMap{h},'PSC'))
                    vecCompPSC{k}=load(num2str(namesMap{h}));
                elseif sum(strfind(namesMap{h},['LSC_',expDate{k}])) && sum(strfind(namesMap{h},'LSC'))
                    vecCompLSC{k}=load(num2str(namesMap{h}));
                end
            end
        end
        
    case 0
        namesASC=getFileList(filedir,['vecComp_',num2str(discPoints),'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_ASC'],0,'anywhere');
        for k=1:size(namesASC,2)
            ind=strfind(namesASC{k},'_');
            expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
            
            namesMap=getFileList(filedir,['vecComp_',num2str(discPoints),'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f')],0,'anywhere');
            for h=1:size(namesMap,2)
                if sum(strfind(namesMap{h},['ASC_',expDate{k}])) && sum(strfind(namesMap{h},'ASC'))
                    vecCompASC{k}=load(num2str(namesMap{h}));
                elseif sum(strfind(namesMap{h},['PSC_',expDate{k}])) && sum(strfind(namesMap{h},'PSC'))
                    vecCompPSC{k}=load(num2str(namesMap{h}));
                elseif sum(strfind(namesMap{h},['LSC_',expDate{k}])) && sum(strfind(namesMap{h},'LSC'))
                    vecCompLSC{k}=load(num2str(namesMap{h}));
                end
            end
        end
end



for k=1:size(vecCompASC,2)    % number of retinas
    j=0;    % the index in the vecCompAll variables
    
    for j=1:size(vecCompASC{k}.isCART,2)
        if (sum(logical(vecCompASC{k}.ONorONOFForOFF(j)==ONorONOFForOFFstate)) && sum(logical(vecCompASC{k}.isCART(j)==CARTstate))...
                && sum(logical(vecCompASC{k}.isRBPMS(j)==RBPMSstate)) && sum(logical(vecCompASC{k}.isRetro(j)==RETROstate)))
            if sum(RETROstate~=1)
                i=i+1;
                pooledmap.XYZvecComp.X=[pooledmap.XYZvecComp.X;vecCompASC{k}.XYZvecComp.X(j)];
                pooledmap.XYZvecComp.Y=[pooledmap.XYZvecComp.Y;vecCompASC{k}.XYZvecComp.Y(j)];
                pooledmap.XYZvecComp.Z=[pooledmap.XYZvecComp.Z;vecCompASC{k}.XYZvecComp.Z(j)];
                pooledmap.XYZvecComp.VecX=[pooledmap.XYZvecComp.VecX;vecCompASC{k}.XYZvecComp.VecX(j)];
                pooledmap.XYZvecComp.VecY=[pooledmap.XYZvecComp.VecY;vecCompASC{k}.XYZvecComp.VecY(j)];
                pooledmap.XYZvecComp.VecZ=[pooledmap.XYZvecComp.VecZ;vecCompASC{k}.XYZvecComp.VecZ(j)];
                pooledmap.XYZvecComp.VecXcompASC=[pooledmap.XYZvecComp.VecXcompASC;vecCompASC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompASC=[pooledmap.XYZvecComp.VecYcompASC;vecCompASC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompASC=[pooledmap.XYZvecComp.VecZcompASC;vecCompASC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZASC=[pooledmap.XYZvecComp.angleXYZASC;vecCompASC{k}.XYZvecComp.angleXYZ(j)];
                pooledmap.XYZvecComp.VecXcompPSC=[pooledmap.XYZvecComp.VecXcompPSC;vecCompPSC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompPSC=[pooledmap.XYZvecComp.VecYcompPSC;vecCompPSC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompPSC=[pooledmap.XYZvecComp.VecZcompPSC;vecCompPSC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZPSC=[pooledmap.XYZvecComp.angleXYZPSC;vecCompPSC{k}.XYZvecComp.angleXYZ(j)];
                pooledmap.XYZvecComp.VecXcompLSC=[pooledmap.XYZvecComp.VecXcompLSC;vecCompLSC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompLSC=[pooledmap.XYZvecComp.VecYcompLSC;vecCompLSC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompLSC=[pooledmap.XYZvecComp.VecZcompLSC;vecCompLSC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZLSC=[pooledmap.XYZvecComp.angleXYZLSC;vecCompLSC{k}.XYZvecComp.angleXYZ(j)];
                
                pooledmap.UVvecComp.U=[pooledmap.UVvecComp.U;vecCompASC{k}.UVvecComp.U(j)];
                pooledmap.UVvecComp.V=[pooledmap.UVvecComp.V;vecCompASC{k}.UVvecComp.V(j)];
                pooledmap.UVvecComp.VecU=[pooledmap.UVvecComp.VecU;vecCompASC{k}.UVvecComp.VecU(j)];
                pooledmap.UVvecComp.VecV=[pooledmap.UVvecComp.VecV;vecCompASC{k}.UVvecComp.VecV(j)];
                pooledmap.UVvecComp.VecUcompASC=[pooledmap.UVvecComp.VecUcompASC;vecCompASC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompASC=[pooledmap.UVvecComp.VecVcompASC;vecCompASC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVASC=[pooledmap.UVvecComp.angleUVASC;vecCompASC{k}.UVvecComp.angleUV(j)];
                pooledmap.UVvecComp.VecUcompPSC=[pooledmap.UVvecComp.VecUcompPSC;vecCompPSC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompPSC=[pooledmap.UVvecComp.VecVcompPSC;vecCompPSC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVPSC=[pooledmap.UVvecComp.angleUVPSC;vecCompPSC{k}.UVvecComp.angleUV(j)];
                pooledmap.UVvecComp.VecUcompLSC=[pooledmap.UVvecComp.VecUcompLSC;vecCompLSC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompLSC=[pooledmap.UVvecComp.VecVcompLSC;vecCompLSC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVLSC=[pooledmap.UVvecComp.angleUVLSC;vecCompLSC{k}.UVvecComp.angleUV(j)];
                
                pooledmap.UVStvecComp.USt=[pooledmap.UVStvecComp.USt;vecCompASC{k}.UVStvecComp.USt(j)];
                pooledmap.UVStvecComp.VSt=[pooledmap.UVStvecComp.VSt;vecCompASC{k}.UVStvecComp.VSt(j)];
                pooledmap.UVStvecComp.VecUSt=[pooledmap.UVStvecComp.VecUSt;vecCompASC{k}.UVStvecComp.VecUSt(j)];
                pooledmap.UVStvecComp.VecVSt=[pooledmap.UVStvecComp.VecVSt;vecCompASC{k}.UVStvecComp.VecVSt(j)];
                pooledmap.UVStvecComp.VecUStcompASC=[pooledmap.UVStvecComp.VecUStcompASC;vecCompASC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompASC=[pooledmap.UVStvecComp.VecVStcompASC;vecCompASC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStASC=[pooledmap.UVStvecComp.angleUVStASC;vecCompASC{k}.UVStvecComp.angleUVSt(j)];
                pooledmap.UVStvecComp.VecUStcompPSC=[pooledmap.UVStvecComp.VecUStcompPSC;vecCompPSC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompPSC=[pooledmap.UVStvecComp.VecVStcompPSC;vecCompPSC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStPSC=[pooledmap.UVStvecComp.angleUVStPSC;vecCompPSC{k}.UVStvecComp.angleUVSt(j)];
                pooledmap.UVStvecComp.VecUStcompLSC=[pooledmap.UVStvecComp.VecUStcompLSC;vecCompLSC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompLSC=[pooledmap.UVStvecComp.VecVStcompLSC;vecCompLSC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStLSC=[pooledmap.UVStvecComp.angleUVStLSC;vecCompLSC{k}.UVStvecComp.angleUVSt(j)];
                
                pooledmap.cellID(i,1)=vecCompASC{k}.cellID(j);
                pooledmap.injectionSite(i,1)=vecCompASC{k}.injectionSite(j);
                pooledmap.ONorONOFForOFF(i,1)=vecCompASC{k}.ONorONOFForOFF(j);
                pooledmap.isCART(i,1)=vecCompASC{k}.isCART(j);
                pooledmap.isRBPMS(i,1)=vecCompASC{k}.isRBPMS(j);
                pooledmap.isRetro(i,1)=vecCompASC{k}.isRetro(j);
                pooledmap.pResponse{i,1}=vecCompASC{k}.pResponse{j};
                pooledmap.OSI(i,1)=vecCompASC{k}.OSI(j);
                
           elseif RETROstate==1 && iscell(vecCompASC{k}.injectionSite(j)) && min(ismember(cell2mat(vecCompASC{k}.injectionSite(j)),INJECTstate))
                i=i+1;
                pooledmap.XYZvecComp.X=[pooledmap.XYZvecComp.X;vecCompASC{k}.XYZvecComp.X(j)];
                pooledmap.XYZvecComp.Y=[pooledmap.XYZvecComp.Y;vecCompASC{k}.XYZvecComp.Y(j)];
                pooledmap.XYZvecComp.Z=[pooledmap.XYZvecComp.Z;vecCompASC{k}.XYZvecComp.Z(j)];
                pooledmap.XYZvecComp.VecX=[pooledmap.XYZvecComp.VecX;vecCompASC{k}.XYZvecComp.VecX(j)];
                pooledmap.XYZvecComp.VecY=[pooledmap.XYZvecComp.VecY;vecCompASC{k}.XYZvecComp.VecY(j)];
                pooledmap.XYZvecComp.VecZ=[pooledmap.XYZvecComp.VecZ;vecCompASC{k}.XYZvecComp.VecZ(j)];
                pooledmap.XYZvecComp.VecXcompASC=[pooledmap.XYZvecComp.VecXcompASC;vecCompASC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompASC=[pooledmap.XYZvecComp.VecYcompASC;vecCompASC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompASC=[pooledmap.XYZvecComp.VecZcompASC;vecCompASC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZASC=[pooledmap.XYZvecComp.angleXYZASC;vecCompASC{k}.XYZvecComp.angleXYZ(j)];
                pooledmap.XYZvecComp.VecXcompPSC=[pooledmap.XYZvecComp.VecXcompPSC;vecCompPSC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompPSC=[pooledmap.XYZvecComp.VecYcompPSC;vecCompPSC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompPSC=[pooledmap.XYZvecComp.VecZcompPSC;vecCompPSC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZPSC=[pooledmap.XYZvecComp.angleXYZPSC;vecCompPSC{k}.XYZvecComp.angleXYZ(j)];
                pooledmap.XYZvecComp.VecXcompLSC=[pooledmap.XYZvecComp.VecXcompLSC;vecCompLSC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompLSC=[pooledmap.XYZvecComp.VecYcompLSC;vecCompLSC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompLSC=[pooledmap.XYZvecComp.VecZcompLSC;vecCompLSC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZLSC=[pooledmap.XYZvecComp.angleXYZLSC;vecCompLSC{k}.XYZvecComp.angleXYZ(j)];
                
                pooledmap.UVvecComp.U=[pooledmap.UVvecComp.U;vecCompASC{k}.UVvecComp.U(j)];
                pooledmap.UVvecComp.V=[pooledmap.UVvecComp.V;vecCompASC{k}.UVvecComp.V(j)];
                pooledmap.UVvecComp.VecU=[pooledmap.UVvecComp.VecU;vecCompASC{k}.UVvecComp.VecU(j)];
                pooledmap.UVvecComp.VecV=[pooledmap.UVvecComp.VecV;vecCompASC{k}.UVvecComp.VecV(j)];
                pooledmap.UVvecComp.VecUcompASC=[pooledmap.UVvecComp.VecUcompASC;vecCompASC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompASC=[pooledmap.UVvecComp.VecVcompASC;vecCompASC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVASC=[pooledmap.UVvecComp.angleUVASC;vecCompASC{k}.UVvecComp.angleUV(j)];
                pooledmap.UVvecComp.VecUcompPSC=[pooledmap.UVvecComp.VecUcompPSC;vecCompPSC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompPSC=[pooledmap.UVvecComp.VecVcompPSC;vecCompPSC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVPSC=[pooledmap.UVvecComp.angleUVPSC;vecCompPSC{k}.UVvecComp.angleUV(j)];
                pooledmap.UVvecComp.VecUcompLSC=[pooledmap.UVvecComp.VecUcompLSC;vecCompLSC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompLSC=[pooledmap.UVvecComp.VecVcompLSC;vecCompLSC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVLSC=[pooledmap.UVvecComp.angleUVLSC;vecCompLSC{k}.UVvecComp.angleUV(j)];
                
                pooledmap.UVStvecComp.USt=[pooledmap.UVStvecComp.USt;vecCompASC{k}.UVStvecComp.USt(j)];
                pooledmap.UVStvecComp.VSt=[pooledmap.UVStvecComp.VSt;vecCompASC{k}.UVStvecComp.VSt(j)];
                pooledmap.UVStvecComp.VecUSt=[pooledmap.UVStvecComp.VecUSt;vecCompASC{k}.UVStvecComp.VecUSt(j)];
                pooledmap.UVStvecComp.VecVSt=[pooledmap.UVStvecComp.VecVSt;vecCompASC{k}.UVStvecComp.VecVSt(j)];
                pooledmap.UVStvecComp.VecUStcompASC=[pooledmap.UVStvecComp.VecUStcompASC;vecCompASC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompASC=[pooledmap.UVStvecComp.VecVStcompASC;vecCompASC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStASC=[pooledmap.UVStvecComp.angleUVStASC;vecCompASC{k}.UVStvecComp.angleUVSt(j)];
                pooledmap.UVStvecComp.VecUStcompPSC=[pooledmap.UVStvecComp.VecUStcompPSC;vecCompPSC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompPSC=[pooledmap.UVStvecComp.VecVStcompPSC;vecCompPSC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStPSC=[pooledmap.UVStvecComp.angleUVStPSC;vecCompPSC{k}.UVStvecComp.angleUVSt(j)];
                pooledmap.UVStvecComp.VecUStcompLSC=[pooledmap.UVStvecComp.VecUStcompLSC;vecCompLSC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompLSC=[pooledmap.UVStvecComp.VecVStcompLSC;vecCompLSC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStLSC=[pooledmap.UVStvecComp.angleUVStLSC;vecCompLSC{k}.UVStvecComp.angleUVSt(j)];
                
                pooledmap.cellID(i,1)=vecCompASC{k}.cellID(j);
                pooledmap.injectionSite(i,1)=vecCompASC{k}.injectionSite(j);
                pooledmap.ONorONOFForOFF(i,1)=vecCompASC{k}.ONorONOFForOFF(j);
                pooledmap.isCART(i,1)=vecCompASC{k}.isCART(j);
                pooledmap.isRBPMS(i,1)=vecCompASC{k}.isRBPMS(j);
                pooledmap.isRetro(i,1)=vecCompASC{k}.isRetro(j);
                pooledmap.pResponse{i,1}=vecCompASC{k}.pResponse{j};
                pooledmap.OSI(i,1)=vecCompASC{k}.OSI(j);
                
elseif RETROstate==1 && ~iscell(vecCompASC{k}.injectionSite(j)) && sum(logical(vecCompASC{k}.injectionSite(j)==INJECTstate))
                i=i+1;
                pooledmap.XYZvecComp.X=[pooledmap.XYZvecComp.X;vecCompASC{k}.XYZvecComp.X(j)];
                pooledmap.XYZvecComp.Y=[pooledmap.XYZvecComp.Y;vecCompASC{k}.XYZvecComp.Y(j)];
                pooledmap.XYZvecComp.Z=[pooledmap.XYZvecComp.Z;vecCompASC{k}.XYZvecComp.Z(j)];
                pooledmap.XYZvecComp.VecX=[pooledmap.XYZvecComp.VecX;vecCompASC{k}.XYZvecComp.VecX(j)];
                pooledmap.XYZvecComp.VecY=[pooledmap.XYZvecComp.VecY;vecCompASC{k}.XYZvecComp.VecY(j)];
                pooledmap.XYZvecComp.VecZ=[pooledmap.XYZvecComp.VecZ;vecCompASC{k}.XYZvecComp.VecZ(j)];
                pooledmap.XYZvecComp.VecXcompASC=[pooledmap.XYZvecComp.VecXcompASC;vecCompASC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompASC=[pooledmap.XYZvecComp.VecYcompASC;vecCompASC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompASC=[pooledmap.XYZvecComp.VecZcompASC;vecCompASC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZASC=[pooledmap.XYZvecComp.angleXYZASC;vecCompASC{k}.XYZvecComp.angleXYZ(j)];
                pooledmap.XYZvecComp.VecXcompPSC=[pooledmap.XYZvecComp.VecXcompPSC;vecCompPSC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompPSC=[pooledmap.XYZvecComp.VecYcompPSC;vecCompPSC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompPSC=[pooledmap.XYZvecComp.VecZcompPSC;vecCompPSC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZPSC=[pooledmap.XYZvecComp.angleXYZPSC;vecCompPSC{k}.XYZvecComp.angleXYZ(j)];
                pooledmap.XYZvecComp.VecXcompLSC=[pooledmap.XYZvecComp.VecXcompLSC;vecCompLSC{k}.XYZvecComp.VecXcomp(j)];
                pooledmap.XYZvecComp.VecYcompLSC=[pooledmap.XYZvecComp.VecYcompLSC;vecCompLSC{k}.XYZvecComp.VecYcomp(j)];
                pooledmap.XYZvecComp.VecZcompLSC=[pooledmap.XYZvecComp.VecZcompLSC;vecCompLSC{k}.XYZvecComp.VecZcomp(j)];
                pooledmap.XYZvecComp.angleXYZLSC=[pooledmap.XYZvecComp.angleXYZLSC;vecCompLSC{k}.XYZvecComp.angleXYZ(j)];
                
                pooledmap.UVvecComp.U=[pooledmap.UVvecComp.U;vecCompASC{k}.UVvecComp.U(j)];
                pooledmap.UVvecComp.V=[pooledmap.UVvecComp.V;vecCompASC{k}.UVvecComp.V(j)];
                pooledmap.UVvecComp.VecU=[pooledmap.UVvecComp.VecU;vecCompASC{k}.UVvecComp.VecU(j)];
                pooledmap.UVvecComp.VecV=[pooledmap.UVvecComp.VecV;vecCompASC{k}.UVvecComp.VecV(j)];
                pooledmap.UVvecComp.VecUcompASC=[pooledmap.UVvecComp.VecUcompASC;vecCompASC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompASC=[pooledmap.UVvecComp.VecVcompASC;vecCompASC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVASC=[pooledmap.UVvecComp.angleUVASC;vecCompASC{k}.UVvecComp.angleUV(j)];
                pooledmap.UVvecComp.VecUcompPSC=[pooledmap.UVvecComp.VecUcompPSC;vecCompPSC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompPSC=[pooledmap.UVvecComp.VecVcompPSC;vecCompPSC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVPSC=[pooledmap.UVvecComp.angleUVPSC;vecCompPSC{k}.UVvecComp.angleUV(j)];
                pooledmap.UVvecComp.VecUcompLSC=[pooledmap.UVvecComp.VecUcompLSC;vecCompLSC{k}.UVvecComp.VecUcomp(j)];
                pooledmap.UVvecComp.VecVcompLSC=[pooledmap.UVvecComp.VecVcompLSC;vecCompLSC{k}.UVvecComp.VecVcomp(j)];
                pooledmap.UVvecComp.angleUVLSC=[pooledmap.UVvecComp.angleUVLSC;vecCompLSC{k}.UVvecComp.angleUV(j)];
                
                pooledmap.UVStvecComp.USt=[pooledmap.UVStvecComp.USt;vecCompASC{k}.UVStvecComp.USt(j)];
                pooledmap.UVStvecComp.VSt=[pooledmap.UVStvecComp.VSt;vecCompASC{k}.UVStvecComp.VSt(j)];
                pooledmap.UVStvecComp.VecUSt=[pooledmap.UVStvecComp.VecUSt;vecCompASC{k}.UVStvecComp.VecUSt(j)];
                pooledmap.UVStvecComp.VecVSt=[pooledmap.UVStvecComp.VecVSt;vecCompASC{k}.UVStvecComp.VecVSt(j)];
                pooledmap.UVStvecComp.VecUStcompASC=[pooledmap.UVStvecComp.VecUStcompASC;vecCompASC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompASC=[pooledmap.UVStvecComp.VecVStcompASC;vecCompASC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStASC=[pooledmap.UVStvecComp.angleUVStASC;vecCompASC{k}.UVStvecComp.angleUVSt(j)];
                pooledmap.UVStvecComp.VecUStcompPSC=[pooledmap.UVStvecComp.VecUStcompPSC;vecCompPSC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompPSC=[pooledmap.UVStvecComp.VecVStcompPSC;vecCompPSC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStPSC=[pooledmap.UVStvecComp.angleUVStPSC;vecCompPSC{k}.UVStvecComp.angleUVSt(j)];
                pooledmap.UVStvecComp.VecUStcompLSC=[pooledmap.UVStvecComp.VecUStcompLSC;vecCompLSC{k}.UVStvecComp.VecUStcomp(j)];
                pooledmap.UVStvecComp.VecVStcompLSC=[pooledmap.UVStvecComp.VecVStcompLSC;vecCompLSC{k}.UVStvecComp.VecVStcomp(j)];
                pooledmap.UVStvecComp.angleUVStLSC=[pooledmap.UVStvecComp.angleUVStLSC;vecCompLSC{k}.UVStvecComp.angleUVSt(j)];
                
                pooledmap.cellID(i,1)=vecCompASC{k}.cellID(j);
                pooledmap.injectionSite(i,1)=vecCompASC{k}.injectionSite(j);
                pooledmap.ONorONOFForOFF(i,1)=vecCompASC{k}.ONorONOFForOFF(j);
                pooledmap.isCART(i,1)=vecCompASC{k}.isCART(j);
                pooledmap.isRBPMS(i,1)=vecCompASC{k}.isRBPMS(j);
                pooledmap.isRetro(i,1)=vecCompASC{k}.isRetro(j);
                pooledmap.pResponse{i,1}=vecCompASC{k}.pResponse{j};
                pooledmap.OSI(i,1)=vecCompASC{k}.OSI(j);
                
            end
        end
    end
end   % for analyzing all the retinas

% figure;
% scale=0.2;
% quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,pooledmap.XYZvecComp.VecX,pooledmap.XYZvecComp.VecY,pooledmap.XYZvecComp.VecZ,scale);
% figure;
% scale=0.2;
% quiver(pooledmap.UVvecComp.U,pooledmap.UVvecComp.V,pooledmap.UVvecComp.VecU,pooledmap.UVvecComp.VecV,scale);


pooledmap.metadata.CARTstate=CARTstate;
pooledmap.metadata.RBPMSstate=RBPMSstate;
pooledmap.metadata.RETROstate=RETROstate;
pooledmap.metadata.INJECTstate=INJECTstate;
pooledmap.metadata.ONorONOFForOFFstate=ONorONOFForOFFstate;
pooledmap.metadata.BGstate=BGstate;
pooledmap.metadata.COLORstate=COLORstate;
pooledmap.metadata.discPoints=discPoints;
pooledmap.metadata.fieldType=fieldType;
pooledmap.metadata.namesMap=namesMap;




%%%%%% plot standard retina surface %%%%%%%%
f1=figure;
StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

n=size(StSec1Data,1);

StRHO1=StSec1Data(1:n,1:n);
StF1=StSec1Data(1:n,n+1:2*n);
StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);

StRHO2=StSec2Data(1:n,1:n);
StF2=StSec2Data(1:n,n+1:2*n);
StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);

StRHO3=StSec3Data(1:n,1:n);
StF3=StSec3Data(1:n,n+1:2*n);
StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);

StRHO4=StSec4Data(1:n,1:n);
StF4=StSec4Data(1:n,n+1:2*n);
StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);

hold on;
s1=surf(StU1,StV1,zeros(size(StU1,1),size(StU1,1)),'FaceColor','k','EdgeColor', 'none');
s2=surf(StU2,StV2,zeros(size(StU2,1),size(StU2,1)),'FaceColor','k','EdgeColor', 'none');
s3=surf(StU3,StV3,zeros(size(StU3,1),size(StU3,1)),'FaceColor','k','EdgeColor', 'none');
s4=surf(StU4,StV4,zeros(size(StU4,1),size(StU4,1)), 'FaceColor','k','EdgeColor', 'none');
alpha(.15);
uistack(s1,'bottom');uistack(s2,'bottom');uistack(s3,'bottom');uistack(s4,'bottom');
hold on;

pooledmap.standardRetina.StU1=StU1;
pooledmap.standardRetina.StU2=StU2;
pooledmap.standardRetina.StU3=StU3;
pooledmap.standardRetina.StU4=StU4;
pooledmap.standardRetina.StV1=StV1;
pooledmap.standardRetina.StV2=StV2;
pooledmap.standardRetina.StV3=StV3;
pooledmap.standardRetina.StV4=StV4;


%%%%%%% plot field on standard retina %%%%%%%% 

if strcmp(fieldType,'ASC')
        normal=[0.592,0.764,0.247];
        graphcolor='b';
elseif strcmp(fieldType,'PSC')
        normal=[0.470,-0.764,0.438];
        graphcolor='m';
elseif strcmp(fieldType,'LSC')
        normal=[0.421,-0.056,-0.901];
        graphcolor='r';
end
try VectorCurveStandard(fieldType,normal,discPoints,BGstate,graphcolor,visAngle); end


if strcmp(fieldType,'all')
    normal=[0.592,0.764,0.247];     %ASC
    try VectorCurveStandard(fieldType,normal,discPoints,BGstate,'b',visAngle); end
    normal=[0.470,-0.764,0.438];        %PSC
    try VectorCurveStandard(fieldType,normal,discPoints,BGstate,'m',visAngle); end
    normal=[0.421,-0.056,-0.901];       %LSC
    try VectorCurveStandard(fieldType,normal,discPoints,BGstate,'r',visAngle); end
end


f2=copyobj(gcf,0);  %duplicate figure f1 for polar plots overlay f1
set(f2,'position',[450 150 950 800]);
f3=copyobj(gcf,0);  %duplicate figure f1 for polar plots overlay f1
set(f3,'position',[450 150 950 800]);


set(0,'CurrentFigure', f1) %make the UVstandard figure active to plot the vectors
scale=0.7;
angleUVSt=[rad2deg(pooledmap.UVStvecComp.angleUVStASC),rad2deg(pooledmap.UVStvecComp.angleUVStPSC),rad2deg(pooledmap.UVStvecComp.angleUVStLSC)];
angleUVStInv=abs(angleUVSt-180);    % calculate the inverse angles
if COLORstate==1
    [Ms,Is]=min(angleUVSt,[],2);      % finds the canal that fit the data best
    quiver(pooledmap.UVStvecComp.USt(Is==1),pooledmap.UVStvecComp.VSt(Is==1),pooledmap.UVStvecComp.VecUSt(Is==1),pooledmap.UVStvecComp.VecVSt(Is==1),scale,'b');
    quiver(pooledmap.UVStvecComp.USt(Is==2),pooledmap.UVStvecComp.VSt(Is==2),pooledmap.UVStvecComp.VecUSt(Is==2),pooledmap.UVStvecComp.VecVSt(Is==2),scale,'m');
    quiver(pooledmap.UVStvecComp.USt(Is==3),pooledmap.UVStvecComp.VSt(Is==3),pooledmap.UVStvecComp.VecUSt(Is==3),pooledmap.UVStvecComp.VecVSt(Is==3),scale,'r');
elseif COLORstate==2 
    angleUVStAll=[rad2deg(pooledmap.UVStvecComp.angleUVStASC),angleUVStInv(:,1),rad2deg(pooledmap.UVStvecComp.angleUVStPSC),...
        angleUVStInv(:,2),rad2deg(pooledmap.UVStvecComp.angleUVStLSC),angleUVStInv(:,3)];
    [Mall,Iall]=min(angleUVStAll,[],2);      % finds the canal that fit either the originl data or the inverse data best
    quiver(pooledmap.UVStvecComp.USt(Iall==1),pooledmap.UVStvecComp.VSt(Iall==1),pooledmap.UVStvecComp.VecUSt(Iall==1),pooledmap.UVStvecComp.VecVSt(Iall==1),scale,'b');
    quiver(pooledmap.UVStvecComp.USt(Iall==2),pooledmap.UVStvecComp.VSt(Iall==2),pooledmap.UVStvecComp.VecUSt(Iall==2),pooledmap.UVStvecComp.VecVSt(Iall==2),scale,'b');
    quiver(pooledmap.UVStvecComp.USt(Iall==3),pooledmap.UVStvecComp.VSt(Iall==3),pooledmap.UVStvecComp.VecUSt(Iall==3),pooledmap.UVStvecComp.VecVSt(Iall==3),scale,'m');
    quiver(pooledmap.UVStvecComp.USt(Iall==4),pooledmap.UVStvecComp.VSt(Iall==4),pooledmap.UVStvecComp.VecUSt(Iall==4),pooledmap.UVStvecComp.VecVSt(Iall==4),scale,'m');
    quiver(pooledmap.UVStvecComp.USt(Iall==5),pooledmap.UVStvecComp.VSt(Iall==5),pooledmap.UVStvecComp.VecUSt(Iall==5),pooledmap.UVStvecComp.VecVSt(Iall==5),scale,'r');
    quiver(pooledmap.UVStvecComp.USt(Iall==6),pooledmap.UVStvecComp.VSt(Iall==6),pooledmap.UVStvecComp.VecUSt(Iall==6),pooledmap.UVStvecComp.VecVSt(Iall==6),scale,'r');
elseif COLORstate==3
    quiver(pooledmap.UVStvecComp.USt,pooledmap.UVStvecComp.VSt,pooledmap.UVStvecComp.VecUSt,pooledmap.UVStvecComp.VecVSt,scale,'k','LineWidth',1);
end
set(f1,'position',[450 150 950 800]);
set(gcf, 'color', [1 1 1]);
shg;
set(findobj(gcf, 'type','axes'), 'Visible','off')

clusteringFlag=0;
if CLUSTERstate==2
pooledmap=addClusterGrid(pooledmap,f1,f2,f3);

elseif CLUSTERstate==1
pause(10);
clusteringFlag=get(handles.doneClusteringButton,'value');
while clusteringFlag==0
    pooledmap=addCluster(pooledmap,f1,f2,f3);
    pause(5);
    clusteringFlag=get(handles.doneClusteringButton,'value');
end
end

filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');

switch ALPHACORRstate
    case 1
        savefig(f1,['pooledMapAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f2,['pooledMapAlphaCorrClusters_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f3,['pooledMapAlphaCorrPolarHist_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        save(['pooledMapAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.mat'],'pooledmap');
    case 0
        saveas(f1,['pooledMap_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f2,['pooledMapClusters_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f3,['pooledMapPolarHist_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        save(['pooledMap_',num2str(filename),'_',num2str(discPoints),'_DSI',num2str(OSIlow, '%0.2f'),'_',num2str(OSIhigh, '%0.2f'),'_created_',currenttime,'.mat'],'pooledmap');
end

set(handles.alphaCorrbutton,'value',0);
set(handles.doneClusteringButton,'value',0);
clusteringFlag=0;
clear global"


