function varargout = pooledDataMappingGUI(varargin)
% POOLEDDATAMAPPINGGUI MATLAB code for pooledDataMappingGUI.fig
%      POOLEDDATAMAPPINGGUI, by itself, creates a new POOLEDDATAMAPPINGGUI or raises the existing
%      singleton*.
%
%      H = POOLEDDATAMAPPINGGUI returns the handle to a new POOLEDDATAMAPPINGGUI or the handle to
%      the existing singleton*.
%
%      POOLEDDATAMAPPINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POOLEDDATAMAPPINGGUI.M with the given input arguments.
%
%      POOLEDDATAMAPPINGGUI('Property','Value',...) creates a new POOLEDDATAMAPPINGGUI or raises the
%      existing singleton*.  Starting from the left, propecrty value pairs are
%      applied to the GUI before pooledDataMappingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pooledDataMappingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pooledDataMappingGUI

% Last Modified by GUIDE v2.5 14-May-2015 10:29:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pooledDataMappingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @pooledDataMappingGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pooledDataMappingGUI is made visible.
function pooledDataMappingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pooledDataMappingGUI (see VARARGIN)

% Choose default command line output for pooledDataMappingGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pooledDataMappingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% global CARTbutton
% global RBPMSbutton
% global discPoints
% 
% CARTbutton=get(handles.CARTbutton,'value');
% RBPMSbutton=get(handles.RBPMSbutton,'value');
% discPoints=str2num(get(handles.discPoints,'String'));



%%%%%% initialize_input_parameters
clear global"

global fieldType
global visAngle
global CARTstate
global RBPMSstate
global ONorONOFForOFFstate
global RETROstate
global INJECTstate
global BGstate
global COLORstate
global CLUSTERstate
global ALPHACORRstate
global DATAMODELstate
global resolution
global polarity2
global DSIlow
global DSIhigh

%initialize handles and display default parameters on user interface
discPoints=40;
visAngle=180;
DSIlow=0.2;
DSIhigh=1;
fieldType='ASC';
resolution=10;
set(handles.ONDSbutton,'value',1);
ONorONOFForOFFstate=1;
set(handles.MTNDbutton,'value',1);
INJECTstate=1;
set(handles.CARTbutton,'value',1);
CARTstate=1;
set(handles.RBPMSbutton,'value',1);
RBPMSstate=1;
set(handles.RETRObutton,'value',1);
RETROstate=1;
set(handles.showVecFieldbutton,'value',1);
BGstate=1;
set(handles.senseColorbutton,'value',1);
COLORstate=1;
set(handles.manualClusterbutton,'value',1);
CLUSTERstate=1;
set(handles.doneClusteringButton,'value',0);
ALPHACORRstate=0;
DATAMODELstate=0;
set(handles.originalPolbutton,'value',1);
polarity2=1;
set(handles.alphaCorrbutton,'value',0);
set(handles.discPoints,'String', num2str(discPoints)); 
set(handles.visAngle,'String', num2str(visAngle));
set(handles.resolution,'String', num2str(resolution));
set(handles.DSIlow,'String', num2str(DSIlow));
set(handles.DSIhigh,'String', num2str(DSIhigh));

% --- Outputs from this function are returned to the command line.
function varargout = pooledDataMappingGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in startAnalysis.
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
global DSIlow
global DSIhigh

discPoints=str2num(get(handles.discPoints,'String'));

% determine the current ONorONOFForOFFstate
ONorONOFForOFFstate=[];
if get(handles.ONDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,1]; end      % ONDS
if get(handles.ONOFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,2]; end   %ONOFFDS
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
        namesASC=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_ASC'],0,'anywhere');
        
        for k=1:size(namesASC,2)
            ind=strfind(namesASC{k},'_');
            expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
            
            namesMap=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f')],0,'anywhere');
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
        namesASC=getFileList(filedir,['vecComp_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_ASC'],0,'anywhere');
        for k=1:size(namesASC,2)
            ind=strfind(namesASC{k},'_');
            expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
            
            namesMap=getFileList(filedir,['vecComp_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f')],0,'anywhere');
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
                pooledmap.DSI(i,1)=vecCompASC{k}.DSI(j);
                
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
                pooledmap.DSI(i,1)=vecCompASC{k}.DSI(j);
                
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
                pooledmap.DSI(i,1)=vecCompASC{k}.DSI(j);
                
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
        savefig(f1,['pooledMapAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f2,['pooledMapAlphaCorrClusters_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f3,['pooledMapAlphaCorrPolarHist_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        save(['pooledMapAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.mat'],'pooledmap');
    case 0
        saveas(f1,['pooledMap_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f2,['pooledMapClusters_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        saveas(f3,['pooledMapPolarHist_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        save(['pooledMap_',num2str(filename),'_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.mat'],'pooledmap');
end

set(handles.alphaCorrbutton,'value',0);
set(handles.doneClusteringButton,'value',0);
clusteringFlag=0;
clear global"


function discPoints_Callback(hObject, eventdata, handles)
% hObject    handle to discPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of discPoints as text
%        str2double(get(hObject,'String')) returns contents of discPoints as a double

global discPoints

temp_discPoints=str2num(get(handles.discPoints,'String'));
if temp_discPoints<10 || temp_discPoints>50
    h=msgbox('number of discretization points must be between 10 and 50', 'Error','error');
else
    discPoints=temp_discPoints;
end


% --- Executes during object creation, after setting all properties.
function discPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to discPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ONDSbutton.
function ONDSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ONDSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ONDSbutton

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.ONDSbutton,'value',1);
else
     set(handles.ONDSbutton,'value',0)
end


% --- Executes on button press in ONOFFDSbutton.
function ONOFFDSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ONOFFDSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ONOFFDSbutton

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.ONOFFDSbutton,'value',1);
else
     set(handles.ONOFFDSbutton,'value',0)
end



% --- Executes on button press in OFFDSbutton.
function OFFDSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OFFDSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OFFDSbutton

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.OFFDSbutton,'value',1);
else
     set(handles.OFFDSbutton,'value',0)
end


% --- Executes on button press in ONtDSbutton.
function ONtDSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ONtDSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ONtDSbutton

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.ONtDSbutton,'value',1);
else
     set(handles.ONtDSbutton,'value',0)
end





% --- Executes on button press in MTNDbutton.
function MTNDbutton_Callback(hObject, eventdata, handles)
% hObject    handle to MTNDbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MTNDbutton

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.MTNDbutton,'value',1);
else
     set(handles.MTNDbutton,'value',0);
end


% --- Executes on button press in SFbutton.
function SFbutton_Callback(hObject, eventdata, handles)
% hObject    handle to SFbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SFbutton

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.SFbutton,'value',1);
else
     set(handles.SFbutton,'value',0);
end


% --- Executes on button press in NOTbutton.
function NOTbutton_Callback(hObject, eventdata, handles)
% hObject    handle to NOTbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NOTbutton

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.NOTbutton,'value',1);
else
     set(handles.NOTbutton,'value',0);
end


% --- Executes on button press in SCbutton.
function SCbutton_Callback(hObject, eventdata, handles)
% hObject    handle to SCbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SCbutton


if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.SCbutton,'value',1);
else
     set(handles.SCbutton,'value',0);
end


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

global CARTstate

CARTimmunoButton=get(eventdata.NewValue,'tag');
switch CARTimmunoButton
     case 'CARTbutton'
         CARTstate=1;
     case 'nonCARTbutton'
         CARTstate=2;
    case 'allCARTbutton'
         CARTstate=[0,1,2];
end


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global RBPMSstate

RBPMSimmunoButton=get(eventdata.NewValue,'tag');
switch RBPMSimmunoButton
     case 'RBPMSbutton'
         RBPMSstate=1;
     case 'nonRBPMSbutton'
         RBPMSstate=2;
    case 'allRBPMSbutton'
         RBPMSstate=[0,1,2];
end


% --- Executes on button press in doneClusteringButton.
function doneClusteringButton_Callback(hObject, eventdata, handles)
% hObject    handle to doneClusteringButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.doneClusteringButton,'value',1);
else
     set(handles.doneClusteringButton,'value',0);
end


% --- Executes when selected object is changed in uipanel5.
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel5 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

global fieldType

fieldTypeButton=get(eventdata.NewValue,'tag');
switch fieldTypeButton
     case 'ASCbutton'
         fieldType='ASC';
     case 'PSCbutton'
         fieldType='PSC';
     case 'LSCbutton'
         fieldType='LSC';
     case 'allCanalsbutton'
         fieldType='all';
     case 'noFieldbutton'
         fieldType=[];
end


% --- Executes when selected object is changed in uipanel9.
function uipanel9_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel9 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global RETROstate

retroButton=get(eventdata.NewValue,'tag');
switch retroButton
     case 'RETRObutton'
         RETROstate=1;
     case 'nonRETRObutton'
         RETROstate=2;
    case 'allRETRObutton'
         RETROstate=[0,1,2];
end


% --- Executes when selected object is changed in uipanel10.
function uipanel10_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel10 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global BGstate

backgroundButton=get(eventdata.NewValue,'tag');
switch backgroundButton
     case 'showVecFieldbutton'
         BGstate=1;
     case 'showCurveFieldbutton'
         BGstate=2;
    case 'showVecCurveFieldbutton'
         BGstate=3;
end


% --- Executes when selected object is changed in uipanel11.
function uipanel11_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel11 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global COLORstate

colorButton=get(eventdata.NewValue,'tag');
switch colorButton
     case 'senseColorbutton'
         COLORstate=1;
     case 'senseAntisenseColorbutton'
         COLORstate=2;
     case 'allBlackbutton'
         COLORstate=3;
end


% --- Executes when selected object is changed in uipanel12.
function uipanel12_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel12 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global CLUSTERstate

clusterButton=get(eventdata.NewValue,'tag');
switch clusterButton
     case 'manualClusterbutton'
         CLUSTERstate=1;
     case 'gridClusterbutton'
         CLUSTERstate=2;
end



function visAngle_Callback(hObject, eventdata, handles)
% hObject    handle to visAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of visAngle as text
%        str2double(get(hObject,'String')) returns contents of visAngle as a double
global visAngle

temp_visAngle=str2num(get(handles.visAngle,'String'));
if temp_visAngle<90 || temp_visAngle>360
    h=msgbox('number of discretization points must be between 90 and 360', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    visAngle=temp_visAngle;
end


% --- Executes during object creation, after setting all properties.
function visAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alphaCorrbutton.
function alphaCorrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to alphaCorrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alphaCorrbutton

global ALPHACORRstate

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.alphaCorrbutton,'value',1);
     ALPHACORRstate=1;
else
     set(handles.alphaCorrbutton,'value',0);
     ALPHACORRstate=0;
end


% --- Executes on button press in allVecFieldsbutton.
function allVecFieldsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to allVecFieldsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end

histType=5;
calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DATAMODELstate
global CLUSTNstate

contourClustButton=get(eventdata.NewValue,'tag');
switch contourClustButton
     case 'pooledmapbutton'
         DATAMODELstate=0;
         CLUSTNstate=0;
         set(handles.clusterNumber, 'enable', 'off')
     case 'clusterbutton'
         DATAMODELstate=1;
         set(handles.clusterNumber, 'enable', 'on')
end        
        



function clusterNumber_Callback(hObject, eventdata, handles)
% hObject    handle to clusterNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clusterNumber as text
%        str2double(get(hObject,'String')) returns contents of clusterNumber as a double
global CLUSTNstate

clustn = str2double(get(hObject,'string'));
if isnan(clustn)
  errordlg('You must enter a numeric value','Invalid Input','modal')
  uicontrol(hObject)
  return
else
  display(clustn)
  CLUSTNstate=clustn;
end


% --- Executes during object creation, after setting all properties.
function clusterNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fitContourbutton.
function fitContourbutton_Callback(hObject, eventdata, handles)
% hObject    handle to fitContourbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end

switch polarity2
    case 1
        histType=1;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        histType=2;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        histType=3;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        fitSurf2(polarity2);
        
    case 2
        
        histType=1;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        histType=2;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        histType=3;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        histType=-1;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        histType=-2;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        histType=-3;
        calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);
        
        fitSurf2(polarity2,resolution);
end


function resolution_Callback(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resolution as text
%        str2double(get(hObject,'String')) returns contents of resolution as a double

global resolution

resolution = str2double(get(hObject,'string'));
if isnan(resolution)
  errordlg('You must enter a numeric value','Invalid Input','modal')
  uicontrol(hObject)
  return
else
  display(resolution)
end


% --- Executes during object creation, after setting all properties.
function resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in theorContourbutton.
function theorContourbutton_Callback(hObject, eventdata, handles)
% hObject    handle to theorContourbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end


histType=4;
calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution);


% --- Executes when selected object is changed in uibuttongroup2.
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup2 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global polarity2

polButton=get(eventdata.NewValue,'tag');
switch polButton
     case 'originalPolbutton'
         polarity2=1;
         set(handles.clusterNumber, 'enable', 'off')
     case 'bothPolbutton'
         polarity2=2;
         set(handles.clusterNumber, 'enable', 'on')
end  



function DSIlow_Callback(hObject, eventdata, handles)
% hObject    handle to DSIlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DSIlow as text
%        str2double(get(hObject,'String')) returns contents of DSIlow as a double

global DSIlow

temp_DSIlow=str2num(get(handles.DSIlow,'String'));
if temp_DSIlow<0 || temp_DSIlow>1
    h=msgbox('DSI must be between 0 and 1', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    DSIlow=temp_DSIlow;
end



% --- Executes during object creation, after setting all properties.
function DSIlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DSIlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DSIhigh_Callback(hObject, eventdata, handles)
% hObject    handle to DSIhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DSIhigh as text
%        str2double(get(hObject,'String')) returns contents of DSIhigh as a double

global DSIhigh

temp_DSIhigh=str2num(get(handles.DSIhigh,'String'));
if temp_DSIhigh<0 || temp_DSIhigh>1
    h=msgbox('DSI must be between 0 and 1', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    DSIhigh=temp_DSIhigh;
end


% --- Executes during object creation, after setting all properties.
function DSIhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DSIhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ONsusOFFDSbutton.
function ONsusOFFDSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ONsusOFFDSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ONsusOFFDSbutton
if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.ONsusOFFDSbutton,'value',1);
else
     set(handles.ONsusOFFDSbutton,'value',0)
end


% --- Executes on button press in poolAllCellsbutton.
function poolAllCellsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to poolAllCellsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHACORRstate

discPoints=str2num(get(handles.discPoints,'String'));

%initializing variables for (X,Y,Z) coordinates
allCells.XYZvecComp.X=[];
allCells.XYZvecComp.Y=[];
allCells.XYZvecComp.Z=[];
allCells.XYZvecComp.VecX=[];
allCells.XYZvecComp.VecY=[];
allCells.XYZvecComp.VecZ=[];
allCells.XYZvecComp.VecXcompASC=[];
allCells.XYZvecComp.VecYcompASC=[];
allCells.XYZvecComp.VecZcompASC=[];
allCells.XYZvecComp.angleXYZASC=[];
allCells.XYZvecComp.VecXcompPSC=[];
allCells.XYZvecComp.VecYcompPSC=[];
allCells.XYZvecComp.VecZcompPSC=[];
allCells.XYZvecComp.angleXYZPSC=[];
allCells.XYZvecComp.VecXcompLSC=[];
allCells.XYZvecComp.VecYcompLSC=[];
allCells.XYZvecComp.VecZcompLSC=[];
allCells.XYZvecComp.angleXYZLSC=[];
%initializing variables for UV coordinates
allCells.UVvecComp.U=[];
allCells.UVvecComp.V=[];
allCells.UVvecComp.VecU=[];
allCells.UVvecComp.VecV=[];
allCells.UVvecComp.VecUcompASC=[];
allCells.UVvecComp.VecVcompASC=[];
allCells.UVvecComp.angleUVASC=[];
allCells.UVvecComp.VecUcompPSC=[];
allCells.UVvecComp.VecVcompPSC=[];
allCells.UVvecComp.angleUVPSC=[];
allCells.UVvecComp.VecUcompLSC=[];
allCells.UVvecComp.VecVcompLSC=[];
allCells.UVvecComp.angleUVLSC=[];
%initializing variables for UVstandard coordinates
allCells.UVStvecComp.USt=[];
allCells.UVStvecComp.VSt=[];
allCells.UVStvecComp.VecUSt=[];
allCells.UVStvecComp.VecVSt=[];
allCells.UVStvecComp.VecUStcompASC=[];
allCells.UVStvecComp.VecVStcompASC=[];
allCells.UVStvecComp.angleUVStASC=[];
allCells.UVStvecComp.VecUStcompPSC=[];
allCells.UVStvecComp.VecVStcompPSC=[];
allCells.UVStvecComp.angleUVStPSC=[];
allCells.UVStvecComp.VecUStcompLSC=[];
allCells.UVStvecComp.VecVStcompLSC=[];
allCells.UVStvecComp.angleUVStLSC=[];

i=0;        % the index in the allCells variables

filedir=cd;

switch ALPHACORRstate
    case 1
        namesASC=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI0.17_1.00','_ASC'],0,'anywhere');
        
        for k=1:size(namesASC,2)
            ind=strfind(namesASC{k},'_');
            expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
            
            namesMap=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI0.17_1.00'],0,'anywhere');
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
        namesASC=getFileList(filedir,['vecComp_',num2str(discPoints),'_DSI0.17_1.00','_ASC'],0,'anywhere');
        for k=1:size(namesASC,2)
            ind=strfind(namesASC{k},'_');
            expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
            
            namesMap=getFileList(filedir,['vecComp_',num2str(discPoints),'_DSI0.17_1.00'],0,'anywhere');
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
        
        i=i+1;
        allCells.XYZvecComp.X=[allCells.XYZvecComp.X;vecCompASC{k}.XYZvecComp.X(j)];
        allCells.XYZvecComp.Y=[allCells.XYZvecComp.Y;vecCompASC{k}.XYZvecComp.Y(j)];
        allCells.XYZvecComp.Z=[allCells.XYZvecComp.Z;vecCompASC{k}.XYZvecComp.Z(j)];
        allCells.XYZvecComp.VecX=[allCells.XYZvecComp.VecX;vecCompASC{k}.XYZvecComp.VecX(j)];
        allCells.XYZvecComp.VecY=[allCells.XYZvecComp.VecY;vecCompASC{k}.XYZvecComp.VecY(j)];
        allCells.XYZvecComp.VecZ=[allCells.XYZvecComp.VecZ;vecCompASC{k}.XYZvecComp.VecZ(j)];
        allCells.XYZvecComp.VecXcompASC=[allCells.XYZvecComp.VecXcompASC;vecCompASC{k}.XYZvecComp.VecXcomp(j)];
        allCells.XYZvecComp.VecYcompASC=[allCells.XYZvecComp.VecYcompASC;vecCompASC{k}.XYZvecComp.VecYcomp(j)];
        allCells.XYZvecComp.VecZcompASC=[allCells.XYZvecComp.VecZcompASC;vecCompASC{k}.XYZvecComp.VecZcomp(j)];
        allCells.XYZvecComp.angleXYZASC=[allCells.XYZvecComp.angleXYZASC;vecCompASC{k}.XYZvecComp.angleXYZ(j)];
        allCells.XYZvecComp.VecXcompPSC=[allCells.XYZvecComp.VecXcompPSC;vecCompPSC{k}.XYZvecComp.VecXcomp(j)];
        allCells.XYZvecComp.VecYcompPSC=[allCells.XYZvecComp.VecYcompPSC;vecCompPSC{k}.XYZvecComp.VecYcomp(j)];
        allCells.XYZvecComp.VecZcompPSC=[allCells.XYZvecComp.VecZcompPSC;vecCompPSC{k}.XYZvecComp.VecZcomp(j)];
        allCells.XYZvecComp.angleXYZPSC=[allCells.XYZvecComp.angleXYZPSC;vecCompPSC{k}.XYZvecComp.angleXYZ(j)];
        allCells.XYZvecComp.VecXcompLSC=[allCells.XYZvecComp.VecXcompLSC;vecCompLSC{k}.XYZvecComp.VecXcomp(j)];
        allCells.XYZvecComp.VecYcompLSC=[allCells.XYZvecComp.VecYcompLSC;vecCompLSC{k}.XYZvecComp.VecYcomp(j)];
        allCells.XYZvecComp.VecZcompLSC=[allCells.XYZvecComp.VecZcompLSC;vecCompLSC{k}.XYZvecComp.VecZcomp(j)];
        allCells.XYZvecComp.angleXYZLSC=[allCells.XYZvecComp.angleXYZLSC;vecCompLSC{k}.XYZvecComp.angleXYZ(j)];
        
        allCells.UVvecComp.U=[allCells.UVvecComp.U;vecCompASC{k}.UVvecComp.U(j)];
        allCells.UVvecComp.V=[allCells.UVvecComp.V;vecCompASC{k}.UVvecComp.V(j)];
        allCells.UVvecComp.VecU=[allCells.UVvecComp.VecU;vecCompASC{k}.UVvecComp.VecU(j)];
        allCells.UVvecComp.VecV=[allCells.UVvecComp.VecV;vecCompASC{k}.UVvecComp.VecV(j)];
        allCells.UVvecComp.VecUcompASC=[allCells.UVvecComp.VecUcompASC;vecCompASC{k}.UVvecComp.VecUcomp(j)];
        allCells.UVvecComp.VecVcompASC=[allCells.UVvecComp.VecVcompASC;vecCompASC{k}.UVvecComp.VecVcomp(j)];
        allCells.UVvecComp.angleUVASC=[allCells.UVvecComp.angleUVASC;vecCompASC{k}.UVvecComp.angleUV(j)];
        allCells.UVvecComp.VecUcompPSC=[allCells.UVvecComp.VecUcompPSC;vecCompPSC{k}.UVvecComp.VecUcomp(j)];
        allCells.UVvecComp.VecVcompPSC=[allCells.UVvecComp.VecVcompPSC;vecCompPSC{k}.UVvecComp.VecVcomp(j)];
        allCells.UVvecComp.angleUVPSC=[allCells.UVvecComp.angleUVPSC;vecCompPSC{k}.UVvecComp.angleUV(j)];
        allCells.UVvecComp.VecUcompLSC=[allCells.UVvecComp.VecUcompLSC;vecCompLSC{k}.UVvecComp.VecUcomp(j)];
        allCells.UVvecComp.VecVcompLSC=[allCells.UVvecComp.VecVcompLSC;vecCompLSC{k}.UVvecComp.VecVcomp(j)];
        allCells.UVvecComp.angleUVLSC=[allCells.UVvecComp.angleUVLSC;vecCompLSC{k}.UVvecComp.angleUV(j)];
        
        allCells.UVStvecComp.USt=[allCells.UVStvecComp.USt;vecCompASC{k}.UVStvecComp.USt(j)];
        allCells.UVStvecComp.VSt=[allCells.UVStvecComp.VSt;vecCompASC{k}.UVStvecComp.VSt(j)];
        allCells.UVStvecComp.VecUSt=[allCells.UVStvecComp.VecUSt;vecCompASC{k}.UVStvecComp.VecUSt(j)];
        allCells.UVStvecComp.VecVSt=[allCells.UVStvecComp.VecVSt;vecCompASC{k}.UVStvecComp.VecVSt(j)];
        allCells.UVStvecComp.VecUStcompASC=[allCells.UVStvecComp.VecUStcompASC;vecCompASC{k}.UVStvecComp.VecUStcomp(j)];
        allCells.UVStvecComp.VecVStcompASC=[allCells.UVStvecComp.VecVStcompASC;vecCompASC{k}.UVStvecComp.VecVStcomp(j)];
        allCells.UVStvecComp.angleUVStASC=[allCells.UVStvecComp.angleUVStASC;vecCompASC{k}.UVStvecComp.angleUVSt(j)];
        allCells.UVStvecComp.VecUStcompPSC=[allCells.UVStvecComp.VecUStcompPSC;vecCompPSC{k}.UVStvecComp.VecUStcomp(j)];
        allCells.UVStvecComp.VecVStcompPSC=[allCells.UVStvecComp.VecVStcompPSC;vecCompPSC{k}.UVStvecComp.VecVStcomp(j)];
        allCells.UVStvecComp.angleUVStPSC=[allCells.UVStvecComp.angleUVStPSC;vecCompPSC{k}.UVStvecComp.angleUVSt(j)];
        allCells.UVStvecComp.VecUStcompLSC=[allCells.UVStvecComp.VecUStcompLSC;vecCompLSC{k}.UVStvecComp.VecUStcomp(j)];
        allCells.UVStvecComp.VecVStcompLSC=[allCells.UVStvecComp.VecVStcompLSC;vecCompLSC{k}.UVStvecComp.VecVStcomp(j)];
        allCells.UVStvecComp.angleUVStLSC=[allCells.UVStvecComp.angleUVStLSC;vecCompLSC{k}.UVStvecComp.angleUVSt(j)];
        
        allCells.cellID(i,1)=vecCompASC{k}.cellID(j);
        allCells.injectionSite(i,1)=vecCompASC{k}.injectionSite(j);
        allCells.ONorONOFForOFF(i,1)=vecCompASC{k}.ONorONOFForOFF(j);
        allCells.isCART(i,1)=vecCompASC{k}.isCART(j);
        allCells.isRBPMS(i,1)=vecCompASC{k}.isRBPMS(j);
        allCells.isRetro(i,1)=vecCompASC{k}.isRetro(j);
        allCells.pResponse{i,1}=vecCompASC{k}.pResponse{j};
        allCells.DSI(i,1)=vecCompASC{k}.DSI(j);
        
    end
end   % for analyzing all the retinas


filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');

switch ALPHACORRstate
    case 1
          save(['allCellsAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_DSI0.17_1.00','_created_',currenttime,'.mat'],'allCells');
    case 0
          save(['allCells_',num2str(filename),'_',num2str(discPoints),'_DSI0.17_1.00','_created_',currenttime,'.mat'],'allCells');
end
