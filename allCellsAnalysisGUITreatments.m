function varargout = allCellsAnalysisGUITreatments(varargin)
% ALLCELLSANALYSISGUITREATMENTS MATLAB code for allCellsAnalysisGUITreatments.fig
%      ALLCELLSANALYSISGUITREATMENTS, by itself, creates a new ALLCELLSANALYSISGUITREATMENTS or raises the existing
%      singleton*.
%
%      H = ALLCELLSANALYSISGUITREATMENTS returns the handle to a new ALLCELLSANALYSISGUITREATMENTS or the handle to
%      the existing singleton*.
%
%      ALLCELLSANALYSISGUITREATMENTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALLCELLSANALYSISGUITREATMENTS.M with the given input arguments.
%
%      ALLCELLSANALYSISGUITREATMENTS('Property','Value',...) creates a new ALLCELLSANALYSISGUITREATMENTS or raises the
%      existing singleton*.  Starting from the left, propecrty value pairs are
%      applied to the GUI before allCellsAnalysisGUITreatments_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to allCellsAnalysisGUITreatments_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help allCellsAnalysisGUITreatments

% Last Modified by GUIDE v2.5 20-Jan-2016 14:48:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @allCellsAnalysisGUITreatments_OpeningFcn, ...
                   'gui_OutputFcn',  @allCellsAnalysisGUITreatments_OutputFcn, ...
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


% --- Executes just before allCellsAnalysisGUITreatments is made visible.
function allCellsAnalysisGUITreatments_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to allCellsAnalysisGUITreatments (see VARARGIN)

% Choose default command line output for allCellsAnalysisGUITreatments
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes allCellsAnalysisGUITreatments wait for user response (see UIRESUME)
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
global transFieldType
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
global TYPECLUSTstate
global TYPEFEATURECLUSTstate
global TYPESUPERCLUSTstate
global USETYPECLUSTstate
global DATAMODELstate
global resolution
global polarity2
global DSIlow
global DSIhigh
global NORMALRESstate
global postPcut
global numIter
global addNoise
global magCut
global alphashift
global ABstate
global locDist
global widthDist
global DISTstate
global mapType
global GSstate


%initialize handles and display default parameters on user interface
discPoints=40;
visAngle=180;
DSIlow=0.3;
DSIhigh=1;
magCut=0;
postPcut=0.5;
alphashift=0;
fieldType='ASC';
transFieldType='nasal';
resolution=5;
numIter=10;
addNoise=0;
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

set(handles.featureNonSupervisedbutton,'value',1);
USETYPECLUSTstate=3;

set(handles.globalSpacebutton,'value',0);
GSstate=0;


ALPHACORRstate=1;
TYPECLUSTstate=0;
TYPEFEATURECLUSTstate=1;
TYPESUPERCLUSTstate=0;
USETYPECLUSTstate=3;
NORMALRESstate=1;

set(handles.pooledmapbutton,'value',1);
set(handles.alphabetaFromFilebutton,'value',1);
set(handles.alphashift, 'enable', 'off')
ABstate=1;

set(handles.applyDistFuncbutton,'value',1);
set(handles.locDist, 'enable', 'on');
set(handles.widthDist, 'enable', 'on');
locDist=pi/4;
widthDist=pi/6;
set(handles.locDist,'String', num2str(locDist));
set(handles.widthDist,'String', num2str(widthDist));
DISTstate=1;

set(handles.rotbutton,'value',1);
mapType=1;


DATAMODELstate=0;
set(handles.bothPolbutton,'value',1);
polarity2=2;
set(handles.alphaCorrbutton,'value',1);
set(handles.typeClustbutton,'value',0);
set(handles.typeFeatureClustbutton,'value',1);
set(handles.typeSuperClustbutton,'value',0);
set(handles.supervisedbutton,'value',0);
set(handles.normalResbutton,'value',1);
set(handles.discPoints,'String', num2str(discPoints)); 
set(handles.visAngle,'String', num2str(visAngle));
set(handles.resolution,'String', num2str(resolution));
set(handles.DSIlow,'String', num2str(DSIlow));
set(handles.DSIhigh,'String', num2str(DSIhigh));
set(handles.magCut,'String', num2str(magCut));
set(handles.postPcut,'String', num2str(postPcut));
set(handles.numIter,'String', num2str(numIter));
set(handles.addNoise,'String', num2str(addNoise));
set(handles.alphashift,'String', num2str(alphashift));


% --- Outputs from this function are returned to the command line.
function varargout = allCellsAnalysisGUITreatments_OutputFcn(hObject, eventdata, handles) 
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
global transFieldType
global visAngle
global ONorONOFForOFFstate
global CARTstate
global RBPMSstate
global INJECTstate
global RETROstate
global BGstate
global COLORstate
global CLUSTERstate
global USETYPECLUSTstate
global ALPHACORRstate
global DSIlow
global DSIhigh
global magCut
global postPcut
global alphashift
global ABstate
global saculaPolarity

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
if get(handles.IFbutton,'value')==1 INJECTstate=[INJECTstate,4]; end        % IF

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
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a allCellsAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a allCells file');
end
load(FileName);

% % determine the source of alpha beta
% if ABstate==1          % get alpha beta from file
% [allPossibleFieldsName,PathName] = uigetfile('*.mat','Select a allPossibleFields file');
% load(allPossibleFieldsName);
% alphatrans=allPossibleFields.alphasing;
% beta=allPossibleFields.betasing;
% elseif ABstate==0          % use fixed alpha beta
%  alphatrans=[0,90,180,270];
%  beta=[90,90,90,90];
% end





% determine the source of alpha beta
        if ABstate==1          % get alpha beta from allPossibleFields file
            [allPossibleFieldsName,PathName] = uigetfile('*.mat','Select a allPossibleFields translation file of real data');
            load(allPossibleFieldsName);
            alphatrans=allPossibleFields.alphasing;
            betatrans=allPossibleFields.betasing;
        elseif ABstate==0          % use fixed alpha beta
            alphatrans=[0,90,180,270];
            betatrans=[90,90,90,90];    
        elseif ABstate==3          % use alpha beta of sacula     
            switch saculaPolarity
                case 1
                    %normal=[-0.180580836088164,-0.331924861656949,0.925859842444695];      % vector of the sacular macula
                    normal=[-0.119,-0.286,0.951];      % vector of the sacular macula Symmetric
                    [alpharot,betarot]=AB(normal);
                    alphaT=rad2deg(alpharot);
                    beta=rad2deg(betarot);
                case -1
                    %normal=[0.180580836088164,0.331924861656949,-0.925859842444695];       % vector of the inverse sacular macula
                    normal=-[-0.119,-0.286,0.951];      % vector of the inverse sacular macula Symmetric
                    [alpharot,betarot]=AB(normal);
                    alphaT=rad2deg(alpharot);
                    beta=rad2deg(betarot);
                    beta=90-abs(90-beta);
                    if alphaT<=180
                        alphaT=180+alphaT;
                    else alphaT>180
                        alphaT=alphaT-180;
                    end
            end
        
            alphatrans=floor(alphaT);
            betatrans=floor(beta);
            
        end





    


    type=(allCells.DSI>=DSIlow) .* (allCells.DSI<=DSIhigh);      % include all DS cells with DSI between DSI low and DSI high in the analysis
    
    if USETYPECLUSTstate==1
        ONorONOFForOFFtmp=allCells.typeClust;
    elseif USETYPECLUSTstate==2
        ONorONOFForOFFtmp=allCells.superTypeClust;
    elseif USETYPECLUSTstate==0
        ONorONOFForOFFtmp=allCells.ONorONOFForOFF;
    elseif USETYPECLUSTstate==3
        ONorONOFForOFFtmp=allCells.featureTypeClust;
    end
    
    j=0;    % the index in the vecCompAll variables
    for j=1:size(allCells.isCART,1)
        if (sum(logical(ONorONOFForOFFtmp(j)==ONorONOFForOFFstate)) && sum(logical(allCells.isCART(j)==CARTstate))...
                && sum(logical(allCells.isRBPMS(j)==RBPMSstate)) && sum(logical(allCells.isRetro(j)==RETROstate)) && logical(type(j)))
            if sum(RETROstate~=1)
                i=i+1;
                pooledmap.XYZvecComp.X=[pooledmap.XYZvecComp.X;allCells.XYZvecComp.X(j)];
                pooledmap.XYZvecComp.Y=[pooledmap.XYZvecComp.Y;allCells.XYZvecComp.Y(j)];
                pooledmap.XYZvecComp.Z=[pooledmap.XYZvecComp.Z;allCells.XYZvecComp.Z(j)];
                pooledmap.XYZvecComp.VecX=[pooledmap.XYZvecComp.VecX;allCells.XYZvecComp.VecX(j)];
                pooledmap.XYZvecComp.VecY=[pooledmap.XYZvecComp.VecY;allCells.XYZvecComp.VecY(j)];
                pooledmap.XYZvecComp.VecZ=[pooledmap.XYZvecComp.VecZ;allCells.XYZvecComp.VecZ(j)];
                pooledmap.XYZvecComp.VecXcompASC=[pooledmap.XYZvecComp.VecXcompASC;allCells.XYZvecComp.VecXcompASC(j)];
                pooledmap.XYZvecComp.VecYcompASC=[pooledmap.XYZvecComp.VecYcompASC;allCells.XYZvecComp.VecYcompASC(j)];
                pooledmap.XYZvecComp.VecZcompASC=[pooledmap.XYZvecComp.VecZcompASC;allCells.XYZvecComp.VecZcompASC(j)];
                pooledmap.XYZvecComp.angleXYZASC=[pooledmap.XYZvecComp.angleXYZASC;allCells.XYZvecComp.angleXYZASC(j)];
                pooledmap.XYZvecComp.VecXcompPSC=[pooledmap.XYZvecComp.VecXcompPSC;allCells.XYZvecComp.VecXcompPSC(j)];
                pooledmap.XYZvecComp.VecYcompPSC=[pooledmap.XYZvecComp.VecYcompPSC;allCells.XYZvecComp.VecYcompPSC(j)];
                pooledmap.XYZvecComp.VecZcompPSC=[pooledmap.XYZvecComp.VecZcompPSC;allCells.XYZvecComp.VecZcompPSC(j)];
                pooledmap.XYZvecComp.angleXYZPSC=[pooledmap.XYZvecComp.angleXYZPSC;allCells.XYZvecComp.angleXYZPSC(j)];
                pooledmap.XYZvecComp.VecXcompLSC=[pooledmap.XYZvecComp.VecXcompLSC;allCells.XYZvecComp.VecXcompLSC(j)];
                pooledmap.XYZvecComp.VecYcompLSC=[pooledmap.XYZvecComp.VecYcompLSC;allCells.XYZvecComp.VecYcompLSC(j)];
                pooledmap.XYZvecComp.VecZcompLSC=[pooledmap.XYZvecComp.VecZcompLSC;allCells.XYZvecComp.VecZcompLSC(j)];
                pooledmap.XYZvecComp.angleXYZLSC=[pooledmap.XYZvecComp.angleXYZLSC;allCells.XYZvecComp.angleXYZLSC(j)];
                
                pooledmap.UVvecComp.U=[pooledmap.UVvecComp.U;allCells.UVvecComp.U(j)];
                pooledmap.UVvecComp.V=[pooledmap.UVvecComp.V;allCells.UVvecComp.V(j)];
                pooledmap.UVvecComp.VecU=[pooledmap.UVvecComp.VecU;allCells.UVvecComp.VecU(j)];
                pooledmap.UVvecComp.VecV=[pooledmap.UVvecComp.VecV;allCells.UVvecComp.VecV(j)];
                pooledmap.UVvecComp.VecUcompASC=[pooledmap.UVvecComp.VecUcompASC;allCells.UVvecComp.VecUcompASC(j)];
                pooledmap.UVvecComp.VecVcompASC=[pooledmap.UVvecComp.VecVcompASC;allCells.UVvecComp.VecVcompASC(j)];
                pooledmap.UVvecComp.angleUVASC=[pooledmap.UVvecComp.angleUVASC;allCells.UVvecComp.angleUVASC(j)];
                pooledmap.UVvecComp.VecUcompPSC=[pooledmap.UVvecComp.VecUcompPSC;allCells.UVvecComp.VecUcompPSC(j)];
                pooledmap.UVvecComp.VecVcompPSC=[pooledmap.UVvecComp.VecVcompPSC;allCells.UVvecComp.VecVcompPSC(j)];
                pooledmap.UVvecComp.angleUVPSC=[pooledmap.UVvecComp.angleUVPSC;allCells.UVvecComp.angleUVPSC(j)];
                pooledmap.UVvecComp.VecUcompLSC=[pooledmap.UVvecComp.VecUcompLSC;allCells.UVvecComp.VecUcompLSC(j)];
                pooledmap.UVvecComp.VecVcompLSC=[pooledmap.UVvecComp.VecVcompLSC;allCells.UVvecComp.VecVcompLSC(j)];
                pooledmap.UVvecComp.angleUVLSC=[pooledmap.UVvecComp.angleUVLSC;allCells.UVvecComp.angleUVLSC(j)];
                
                pooledmap.UVStvecComp.USt=[pooledmap.UVStvecComp.USt;allCells.UVStvecComp.USt(j)];
                pooledmap.UVStvecComp.VSt=[pooledmap.UVStvecComp.VSt;allCells.UVStvecComp.VSt(j)];
                pooledmap.UVStvecComp.VecUSt=[pooledmap.UVStvecComp.VecUSt;allCells.UVStvecComp.VecUSt(j)];
                pooledmap.UVStvecComp.VecVSt=[pooledmap.UVStvecComp.VecVSt;allCells.UVStvecComp.VecVSt(j)];
                pooledmap.UVStvecComp.VecUStcompASC=[pooledmap.UVStvecComp.VecUStcompASC;allCells.UVStvecComp.VecUStcompASC(j)];
                pooledmap.UVStvecComp.VecVStcompASC=[pooledmap.UVStvecComp.VecVStcompASC;allCells.UVStvecComp.VecVStcompASC(j)];
                pooledmap.UVStvecComp.angleUVStASC=[pooledmap.UVStvecComp.angleUVStASC;allCells.UVStvecComp.angleUVStASC(j)];
                pooledmap.UVStvecComp.VecUStcompPSC=[pooledmap.UVStvecComp.VecUStcompPSC;allCells.UVStvecComp.VecUStcompPSC(j)];
                pooledmap.UVStvecComp.VecVStcompPSC=[pooledmap.UVStvecComp.VecVStcompPSC;allCells.UVStvecComp.VecVStcompPSC(j)];
                pooledmap.UVStvecComp.angleUVStPSC=[pooledmap.UVStvecComp.angleUVStPSC;allCells.UVStvecComp.angleUVStPSC(j)];
                pooledmap.UVStvecComp.VecUStcompLSC=[pooledmap.UVStvecComp.VecUStcompLSC;allCells.UVStvecComp.VecUStcompLSC(j)];
                pooledmap.UVStvecComp.VecVStcompLSC=[pooledmap.UVStvecComp.VecVStcompLSC;allCells.UVStvecComp.VecVStcompLSC(j)];
                pooledmap.UVStvecComp.angleUVStLSC=[pooledmap.UVStvecComp.angleUVStLSC;allCells.UVStvecComp.angleUVStLSC(j)];
                
                pooledmap.cellID(i,1)=allCells.cellID(j);
                pooledmap.injectionSite(i,1)=allCells.injectionSite(j);
                pooledmap.ONorONOFForOFF(i,1)=allCells.ONorONOFForOFF(j);
                try pooledmap.typeClust(i,1)=allCells.typeClust(j); end;
                try pooledmap.superTypeClust(i,1)=allCells.superTypeClust(j); end;
                pooledmap.isCART(i,1)=allCells.isCART(j);
                pooledmap.isRBPMS(i,1)=allCells.isRBPMS(j);
                pooledmap.isRetro(i,1)=allCells.isRetro(j);
                pooledmap.pResponse{i,1}=allCells.pResponse{j};
                pooledmap.DSI(i,1)=allCells.DSI(j);
                pooledmap.onset1(i,1)=allCells.onset1(j);
                
           elseif RETROstate==1 && iscell(allCells.injectionSite(j)) && min(ismember(cell2mat(allCells.injectionSite(j)),INJECTstate))
                i=i+1;
                pooledmap.XYZvecComp.X=[pooledmap.XYZvecComp.X;allCells.XYZvecComp.X(j)];
                pooledmap.XYZvecComp.Y=[pooledmap.XYZvecComp.Y;allCells.XYZvecComp.Y(j)];
                pooledmap.XYZvecComp.Z=[pooledmap.XYZvecComp.Z;allCells.XYZvecComp.Z(j)];
                pooledmap.XYZvecComp.VecX=[pooledmap.XYZvecComp.VecX;allCells.XYZvecComp.VecX(j)];
                pooledmap.XYZvecComp.VecY=[pooledmap.XYZvecComp.VecY;allCells.XYZvecComp.VecY(j)];
                pooledmap.XYZvecComp.VecZ=[pooledmap.XYZvecComp.VecZ;allCells.XYZvecComp.VecZ(j)];
                pooledmap.XYZvecComp.VecXcompASC=[pooledmap.XYZvecComp.VecXcompASC;allCells.XYZvecComp.VecXcompASC(j)];
                pooledmap.XYZvecComp.VecYcompASC=[pooledmap.XYZvecComp.VecYcompASC;allCells.XYZvecComp.VecYcompASC(j)];
                pooledmap.XYZvecComp.VecZcompASC=[pooledmap.XYZvecComp.VecZcompASC;allCells.XYZvecComp.VecZcompASC(j)];
                pooledmap.XYZvecComp.angleXYZASC=[pooledmap.XYZvecComp.angleXYZASC;allCells.XYZvecComp.angleXYZASC(j)];
                pooledmap.XYZvecComp.VecXcompPSC=[pooledmap.XYZvecComp.VecXcompPSC;allCells.XYZvecComp.VecXcompPSC(j)];
                pooledmap.XYZvecComp.VecYcompPSC=[pooledmap.XYZvecComp.VecYcompPSC;allCells.XYZvecComp.VecYcompPSC(j)];
                pooledmap.XYZvecComp.VecZcompPSC=[pooledmap.XYZvecComp.VecZcompPSC;allCells.XYZvecComp.VecZcompPSC(j)];
                pooledmap.XYZvecComp.angleXYZPSC=[pooledmap.XYZvecComp.angleXYZPSC;allCells.XYZvecComp.angleXYZPSC(j)];
                pooledmap.XYZvecComp.VecXcompLSC=[pooledmap.XYZvecComp.VecXcompLSC;allCells.XYZvecComp.VecXcompLSC(j)];
                pooledmap.XYZvecComp.VecYcompLSC=[pooledmap.XYZvecComp.VecYcompLSC;allCells.XYZvecComp.VecYcompLSC(j)];
                pooledmap.XYZvecComp.VecZcompLSC=[pooledmap.XYZvecComp.VecZcompLSC;allCells.XYZvecComp.VecZcompLSC(j)];
                pooledmap.XYZvecComp.angleXYZLSC=[pooledmap.XYZvecComp.angleXYZLSC;allCells.XYZvecComp.angleXYZLSC(j)];
                
                pooledmap.UVvecComp.U=[pooledmap.UVvecComp.U;allCells.UVvecComp.U(j)];
                pooledmap.UVvecComp.V=[pooledmap.UVvecComp.V;allCells.UVvecComp.V(j)];
                pooledmap.UVvecComp.VecU=[pooledmap.UVvecComp.VecU;allCells.UVvecComp.VecU(j)];
                pooledmap.UVvecComp.VecV=[pooledmap.UVvecComp.VecV;allCells.UVvecComp.VecV(j)];
                pooledmap.UVvecComp.VecUcompASC=[pooledmap.UVvecComp.VecUcompASC;allCells.UVvecComp.VecUcompASC(j)];
                pooledmap.UVvecComp.VecVcompASC=[pooledmap.UVvecComp.VecVcompASC;allCells.UVvecComp.VecVcompASC(j)];
                pooledmap.UVvecComp.angleUVASC=[pooledmap.UVvecComp.angleUVASC;allCells.UVvecComp.angleUVASC(j)];
                pooledmap.UVvecComp.VecUcompPSC=[pooledmap.UVvecComp.VecUcompPSC;allCells.UVvecComp.VecUcompPSC(j)];
                pooledmap.UVvecComp.VecVcompPSC=[pooledmap.UVvecComp.VecVcompPSC;allCells.UVvecComp.VecVcompPSC(j)];
                pooledmap.UVvecComp.angleUVPSC=[pooledmap.UVvecComp.angleUVPSC;allCells.UVvecComp.angleUVPSC(j)];
                pooledmap.UVvecComp.VecUcompLSC=[pooledmap.UVvecComp.VecUcompLSC;allCells.UVvecComp.VecUcompLSC(j)];
                pooledmap.UVvecComp.VecVcompLSC=[pooledmap.UVvecComp.VecVcompLSC;allCells.UVvecComp.VecVcompLSC(j)];
                pooledmap.UVvecComp.angleUVLSC=[pooledmap.UVvecComp.angleUVLSC;allCells.UVvecComp.angleUVLSC(j)];
                
                pooledmap.UVStvecComp.USt=[pooledmap.UVStvecComp.USt;allCells.UVStvecComp.USt(j)];
                pooledmap.UVStvecComp.VSt=[pooledmap.UVStvecComp.VSt;allCells.UVStvecComp.VSt(j)];
                pooledmap.UVStvecComp.VecUSt=[pooledmap.UVStvecComp.VecUSt;allCells.UVStvecComp.VecUSt(j)];
                pooledmap.UVStvecComp.VecVSt=[pooledmap.UVStvecComp.VecVSt;allCells.UVStvecComp.VecVSt(j)];
                pooledmap.UVStvecComp.VecUStcompASC=[pooledmap.UVStvecComp.VecUStcompASC;allCells.UVStvecComp.VecUStcompASC(j)];
                pooledmap.UVStvecComp.VecVStcompASC=[pooledmap.UVStvecComp.VecVStcompASC;allCells.UVStvecComp.VecVStcompASC(j)];
                pooledmap.UVStvecComp.angleUVStASC=[pooledmap.UVStvecComp.angleUVStASC;allCells.UVStvecComp.angleUVStASC(j)];
                pooledmap.UVStvecComp.VecUStcompPSC=[pooledmap.UVStvecComp.VecUStcompPSC;allCells.UVStvecComp.VecUStcompPSC(j)];
                pooledmap.UVStvecComp.VecVStcompPSC=[pooledmap.UVStvecComp.VecVStcompPSC;allCells.UVStvecComp.VecVStcompPSC(j)];
                pooledmap.UVStvecComp.angleUVStPSC=[pooledmap.UVStvecComp.angleUVStPSC;allCells.UVStvecComp.angleUVStPSC(j)];
                pooledmap.UVStvecComp.VecUStcompLSC=[pooledmap.UVStvecComp.VecUStcompLSC;allCells.UVStvecComp.VecUStcompLSC(j)];
                pooledmap.UVStvecComp.VecVStcompLSC=[pooledmap.UVStvecComp.VecVStcompLSC;allCells.UVStvecComp.VecVStcompLSC(j)];
                pooledmap.UVStvecComp.angleUVStLSC=[pooledmap.UVStvecComp.angleUVStLSC;allCells.UVStvecComp.angleUVStLSC(j)];
                
                pooledmap.cellID(i,1)=allCells.cellID(j);
                pooledmap.injectionSite(i,1)=allCells.injectionSite(j);
                pooledmap.ONorONOFForOFF(i,1)=allCells.ONorONOFForOFF(j);
                try pooledmap.typeClust(i,1)=allCells.typeClust(j); end;
                try pooledmap.superTypeClust(i,1)=allCells.superTypeClust(j); end;
                pooledmap.isCART(i,1)=allCells.isCART(j);
                pooledmap.isRBPMS(i,1)=allCells.isRBPMS(j);
                pooledmap.isRetro(i,1)=allCells.isRetro(j);
                pooledmap.pResponse{i,1}=allCells.pResponse{j};
                pooledmap.DSI(i,1)=allCells.DSI(j);
                pooledmap.onset1(i,1)=allCells.onset1(j);
                
elseif RETROstate==1 && ~iscell(allCells.injectionSite(j)) && sum(logical(allCells.injectionSite(j)==INJECTstate))
                i=i+1;
                pooledmap.XYZvecComp.X=[pooledmap.XYZvecComp.X;allCells.XYZvecComp.X(j)];
                pooledmap.XYZvecComp.Y=[pooledmap.XYZvecComp.Y;allCells.XYZvecComp.Y(j)];
                pooledmap.XYZvecComp.Z=[pooledmap.XYZvecComp.Z;allCells.XYZvecComp.Z(j)];
                pooledmap.XYZvecComp.VecX=[pooledmap.XYZvecComp.VecX;allCells.XYZvecComp.VecX(j)];
                pooledmap.XYZvecComp.VecY=[pooledmap.XYZvecComp.VecY;allCells.XYZvecComp.VecY(j)];
                pooledmap.XYZvecComp.VecZ=[pooledmap.XYZvecComp.VecZ;allCells.XYZvecComp.VecZ(j)];
                pooledmap.XYZvecComp.VecXcompASC=[pooledmap.XYZvecComp.VecXcompASC;allCells.XYZvecComp.VecXcompASC(j)];
                pooledmap.XYZvecComp.VecYcompASC=[pooledmap.XYZvecComp.VecYcompASC;allCells.XYZvecComp.VecYcompASC(j)];
                pooledmap.XYZvecComp.VecZcompASC=[pooledmap.XYZvecComp.VecZcompASC;allCells.XYZvecComp.VecZcompASC(j)];
                pooledmap.XYZvecComp.angleXYZASC=[pooledmap.XYZvecComp.angleXYZASC;allCells.XYZvecComp.angleXYZASC(j)];
                pooledmap.XYZvecComp.VecXcompPSC=[pooledmap.XYZvecComp.VecXcompPSC;allCells.XYZvecComp.VecXcompPSC(j)];
                pooledmap.XYZvecComp.VecYcompPSC=[pooledmap.XYZvecComp.VecYcompPSC;allCells.XYZvecComp.VecYcompPSC(j)];
                pooledmap.XYZvecComp.VecZcompPSC=[pooledmap.XYZvecComp.VecZcompPSC;allCells.XYZvecComp.VecZcompPSC(j)];
                pooledmap.XYZvecComp.angleXYZPSC=[pooledmap.XYZvecComp.angleXYZPSC;allCells.XYZvecComp.angleXYZPSC(j)];
                pooledmap.XYZvecComp.VecXcompLSC=[pooledmap.XYZvecComp.VecXcompLSC;allCells.XYZvecComp.VecXcompLSC(j)];
                pooledmap.XYZvecComp.VecYcompLSC=[pooledmap.XYZvecComp.VecYcompLSC;allCells.XYZvecComp.VecYcompLSC(j)];
                pooledmap.XYZvecComp.VecZcompLSC=[pooledmap.XYZvecComp.VecZcompLSC;allCells.XYZvecComp.VecZcompLSC(j)];
                pooledmap.XYZvecComp.angleXYZLSC=[pooledmap.XYZvecComp.angleXYZLSC;allCells.XYZvecComp.angleXYZLSC(j)];
                
                pooledmap.UVvecComp.U=[pooledmap.UVvecComp.U;allCells.UVvecComp.U(j)];
                pooledmap.UVvecComp.V=[pooledmap.UVvecComp.V;allCells.UVvecComp.V(j)];
                pooledmap.UVvecComp.VecU=[pooledmap.UVvecComp.VecU;allCells.UVvecComp.VecU(j)];
                pooledmap.UVvecComp.VecV=[pooledmap.UVvecComp.VecV;allCells.UVvecComp.VecV(j)];
                pooledmap.UVvecComp.VecUcompASC=[pooledmap.UVvecComp.VecUcompASC;allCells.UVvecComp.VecUcompASC(j)];
                pooledmap.UVvecComp.VecVcompASC=[pooledmap.UVvecComp.VecVcompASC;allCells.UVvecComp.VecVcompASC(j)];
                pooledmap.UVvecComp.angleUVASC=[pooledmap.UVvecComp.angleUVASC;allCells.UVvecComp.angleUVASC(j)];
                pooledmap.UVvecComp.VecUcompPSC=[pooledmap.UVvecComp.VecUcompPSC;allCells.UVvecComp.VecUcompPSC(j)];
                pooledmap.UVvecComp.VecVcompPSC=[pooledmap.UVvecComp.VecVcompPSC;allCells.UVvecComp.VecVcompPSC(j)];
                pooledmap.UVvecComp.angleUVPSC=[pooledmap.UVvecComp.angleUVPSC;allCells.UVvecComp.angleUVPSC(j)];
                pooledmap.UVvecComp.VecUcompLSC=[pooledmap.UVvecComp.VecUcompLSC;allCells.UVvecComp.VecUcompLSC(j)];
                pooledmap.UVvecComp.VecVcompLSC=[pooledmap.UVvecComp.VecVcompLSC;allCells.UVvecComp.VecVcompLSC(j)];
                pooledmap.UVvecComp.angleUVLSC=[pooledmap.UVvecComp.angleUVLSC;allCells.UVvecComp.angleUVLSC(j)];
                
                pooledmap.UVStvecComp.USt=[pooledmap.UVStvecComp.USt;allCells.UVStvecComp.USt(j)];
                pooledmap.UVStvecComp.VSt=[pooledmap.UVStvecComp.VSt;allCells.UVStvecComp.VSt(j)];
                pooledmap.UVStvecComp.VecUSt=[pooledmap.UVStvecComp.VecUSt;allCells.UVStvecComp.VecUSt(j)];
                pooledmap.UVStvecComp.VecVSt=[pooledmap.UVStvecComp.VecVSt;allCells.UVStvecComp.VecVSt(j)];
                pooledmap.UVStvecComp.VecUStcompASC=[pooledmap.UVStvecComp.VecUStcompASC;allCells.UVStvecComp.VecUStcompASC(j)];
                pooledmap.UVStvecComp.VecVStcompASC=[pooledmap.UVStvecComp.VecVStcompASC;allCells.UVStvecComp.VecVStcompASC(j)];
                pooledmap.UVStvecComp.angleUVStASC=[pooledmap.UVStvecComp.angleUVStASC;allCells.UVStvecComp.angleUVStASC(j)];
                pooledmap.UVStvecComp.VecUStcompPSC=[pooledmap.UVStvecComp.VecUStcompPSC;allCells.UVStvecComp.VecUStcompPSC(j)];
                pooledmap.UVStvecComp.VecVStcompPSC=[pooledmap.UVStvecComp.VecVStcompPSC;allCells.UVStvecComp.VecVStcompPSC(j)];
                pooledmap.UVStvecComp.angleUVStPSC=[pooledmap.UVStvecComp.angleUVStPSC;allCells.UVStvecComp.angleUVStPSC(j)];
                pooledmap.UVStvecComp.VecUStcompLSC=[pooledmap.UVStvecComp.VecUStcompLSC;allCells.UVStvecComp.VecUStcompLSC(j)];
                pooledmap.UVStvecComp.VecVStcompLSC=[pooledmap.UVStvecComp.VecVStcompLSC;allCells.UVStvecComp.VecVStcompLSC(j)];
                pooledmap.UVStvecComp.angleUVStLSC=[pooledmap.UVStvecComp.angleUVStLSC;allCells.UVStvecComp.angleUVStLSC(j)];
                
                pooledmap.cellID(i,1)=allCells.cellID(j);
                pooledmap.injectionSite(i,1)=allCells.injectionSite(j);
                pooledmap.ONorONOFForOFF(i,1)=allCells.ONorONOFForOFF(j);
                try pooledmap.typeClust(i,1)=allCells.typeClust(j); end;
                try pooledmap.superTypeClust(i,1)=allCells.superTypeClust(j); end;
                pooledmap.isCART(i,1)=allCells.isCART(j);
                pooledmap.isRBPMS(i,1)=allCells.isRBPMS(j);
                pooledmap.isRetro(i,1)=allCells.isRetro(j);
                pooledmap.pResponse{i,1}=allCells.pResponse{j};
                pooledmap.DSI(i,1)=allCells.DSI(j);
                pooledmap.onset1(i,1)=allCells.onset1(j);
                
            end
        end
    end

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
pooledmap.metadata.USETYPECLUSTstate=USETYPECLUSTstate;
pooledmap.metadata.BGstate=BGstate;
pooledmap.metadata.COLORstate=COLORstate;
pooledmap.metadata.discPoints=discPoints;
pooledmap.metadata.fieldType=fieldType;
pooledmap.metadata.transFieldType=transFieldType;
pooledmap.metadata.namesMap=FileName;
pooledmap.metadata.alphashift=alphashift;

%% plot standard retina surface %%%%%%%%
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

%% plot rotation field on standard retina %%%%%%%% 
switch ABstate
    case 2
alphaTrans=[];
betaTrans=[];
if strcmp(fieldType,'ASC')
        %normal=[0.596,0.766,0.241];
        normal=[0.604,0.758,0.243]; % ASC symmetric
        graphcolor='b';
elseif strcmp(fieldType,'PSC')
        %normal=[0.471,-0.766,0.438];
        normal=[0.447,-0.76,0.469];    % PSC symmetric
        graphcolor='m';
elseif strcmp(fieldType,'LSC')
        %normal=[0.421,-0.060,-0.905];
        normal=[0.449,-0.085,-0.886];   % LSC symmetric
        graphcolor='r';
end

try VectorCurveStandard(fieldType,normal,alphaTrans,betaTrans,discPoints,BGstate,graphcolor,visAngle); end


if strcmp(fieldType,'all')
    %normal=[0.596,0.766,0.241];         %ASC
    normal=[0.604,0.758,0.243]; %ASC symmetric
    try VectorCurveStandard(fieldType,normal,alphaTrans,betaTrans,discPoints,BGstate,'b',visAngle); end
    %normal=[0.471,-0.766,0.438];        %PSC
    normal=[0.447,-0.76,0.469];    % PSC symmetric
    try VectorCurveStandard(fieldType,normal,alphaTrans,betaTrans,discPoints,BGstate,'m',visAngle); end
    %normal=[0.421,-0.060,-0.905];       %LSC
    normal=[0.449,-0.085,-0.886];   % LSC symmetric
    try VectorCurveStandard(fieldType,normal,alphaTrans,betaTrans,discPoints,BGstate,'r',visAngle); end
end

    case 1  % get alpha beta from file
    normal=[];    
        if strcmp(fieldType,'ASC')
            alphaTrans=alphatrans(1);
            betaTrans=betatrans(1);
            graphcolor='c';
        elseif strcmp(fieldType,'iLSC')
            alphaTrans=alphatrans(2);
            betaTrans=betatrans(2);
            graphcolor='c';
        elseif strcmp(fieldType,'iASC')
            alphaTrans=alphatrans(3);
            betaTrans=betatrans(3);
            graphcolor='c';
        elseif strcmp(fieldType,'LSC')
            alphaTrans=alphatrans(4);
            betaTrans=betatrans(4);
            graphcolor='c';
        end

try VectorCurveStandard(fieldType,normal,alphaTrans,betaTrans,discPoints,BGstate,graphcolor,visAngle); end
if strcmp(fieldType,'all')
%         VectorCurveStandardTrans(fieldType,alphabar, beta,discPoints,plotType,graphcolor,VizAngle)
    try VectorCurveStandard(fieldType,normal,alphatrans(1)+alphashift,betatrans(1),discPoints,BGstate,'b',visAngle); end       %nasal     
    try VectorCurveStandard(fieldType,normal,alphatrans(2)+alphashift,betatrans(2),discPoints,BGstate,'g',visAngle); end       %dorsal  
    try VectorCurveStandard(fieldType,normal,alphatrans(3)+alphashift,betatrans(3),discPoints,BGstate,'r',visAngle); end       %temporal   
    try VectorCurveStandard(fieldType,normal,alphatrans(4)+alphashift,betatrans(4),discPoints,BGstate,'k',visAngle); end       %ventral
end

end

%% plot translation field on standard retina %%%%%%%% 

% create spherical retina figure
f4=figure;
f4=RetinaPlot(f4);

alphaTrans=[];
betaTrans=[];
if strcmp(transFieldType,'nasal')
        alphaTrans=alphatrans(1)+alphashift;
        betaTrans=betatrans(1);
        graphcolor='c';
elseif strcmp(transFieldType,'dorsal')
        alphaTrans=alphatrans(2)+alphashift;
        betaTrans=betatrans(2);
        graphcolor='c';
elseif strcmp(transFieldType,'temporal')
        alphaTrans=alphatrans(3)+alphashift;
        betaTrans=betatrans(3);
        graphcolor='c';
elseif strcmp(transFieldType,'ventral')
        alphaTrans=alphatrans(4)+alphashift;
        betaTrans=betatrans(4);
        graphcolor='c';
end

set(0,'CurrentFigure', f1) %make the UVstandard figure active to select a ROI
try VectorCurveStandardTrans(fieldType,alphaTrans,betaTrans,discPoints,BGstate,graphcolor,visAngle); end
set(0,'CurrentFigure', f4) %make the UVstandard figure active to select a ROI
try VectorCurveSphericalTransRetinaFieldLines(alphaTrans,betaTrans,discPoints,visAngle,graphcolor); end
set(0,'CurrentFigure', f1) %make the UVstandard figure active to select a ROI

if strcmp(transFieldType,'all')  
    try VectorCurveStandardTrans(fieldType,alphatrans(1)+alphashift,betatrans(1),discPoints,BGstate,'b',visAngle); end       %nasal
    set(0,'CurrentFigure', f4) 
    try VectorCurveSphericalTransRetinaFieldLines(alphatrans(1),betatrans(1),discPoints,visAngle,'b'); end
    set(0,'CurrentFigure', f1) 
    try VectorCurveStandardTrans(fieldType,alphatrans(2)+alphashift,betatrans(2),discPoints,BGstate,'g',visAngle); end       %dorsal
    set(0,'CurrentFigure', f4) 
    try VectorCurveSphericalTransRetinaFieldLines(alphatrans(2),betatrans(2),discPoints,visAngle,'g'); end
    set(0,'CurrentFigure', f1) 
    try VectorCurveStandardTrans(fieldType,alphatrans(3)+alphashift,betatrans(3),discPoints,BGstate,'r',visAngle); end       %temporal
    set(0,'CurrentFigure', f4) 
    try VectorCurveSphericalTransRetinaFieldLines(alphatrans(3),betatrans(3),discPoints,visAngle,'r'); end
    set(0,'CurrentFigure', f1) 
    try VectorCurveStandardTrans(fieldType,alphatrans(4)+alphashift,betatrans(4),discPoints,BGstate,'m',visAngle); end       %ventral
    set(0,'CurrentFigure', f4) 
    try VectorCurveSphericalTransRetinaFieldLines(alphatrans(4),betatrans(4),discPoints,visAngle,'m'); end
    set(0,'CurrentFigure', f1) 
end

% plots the position of the intersection of the saccule plane with the retina
%plotintersectionLine(discPoints);


%%
f2=copyobj(gcf,0);  %duplicate figure f1 for polar plots overlay f1
set(f2,'position',[450 150 950 800]);
f3=copyobj(gcf,0);  %duplicate figure f1 for polar plots overlay f1
set(f3,'position',[450 150 950 800]);


set(0,'CurrentFigure', f1) %make the UVstandard figure active to plot the vectors
scale=0.4;
angleUVSt=[rad2deg(pooledmap.UVStvecComp.angleUVStASC),rad2deg(pooledmap.UVStvecComp.angleUVStPSC),rad2deg(pooledmap.UVStvecComp.angleUVStLSC)];
angleUVStInv=abs(angleUVSt-180);    % calculate the inverse angles
if COLORstate==0
    angleUVStTrans=[allPossibleFields.anglediff180{allPossibleFields.Irowmax(1),allPossibleFields.Icolmax(1)},...
        allPossibleFields.anglediff180{allPossibleFields.Irowmax(2),allPossibleFields.Icolmax(2)},...
        allPossibleFields.anglediff180{allPossibleFields.Irowmax(3),allPossibleFields.Icolmax(3)},...
        allPossibleFields.anglediff180{allPossibleFields.Irowmax(4),allPossibleFields.Icolmax(4)}];
    
    [Mtrans,Itrans]=min(angleUVStTrans,[],2);      % finds the channel that fits the four translation channels (imported from allPossibleFields file) best
    pooledmap.UVStvecComp.Itrans=Itrans;
    quiver(pooledmap.UVStvecComp.USt(Itrans==1),pooledmap.UVStvecComp.VSt(Itrans==1),pooledmap.UVStvecComp.VecUSt(Itrans==1),pooledmap.UVStvecComp.VecVSt(Itrans==1),scale,'b');
    quiver(pooledmap.UVStvecComp.USt(Itrans==2),pooledmap.UVStvecComp.VSt(Itrans==2),pooledmap.UVStvecComp.VecUSt(Itrans==2),pooledmap.UVStvecComp.VecVSt(Itrans==2),scale,'g');
    quiver(pooledmap.UVStvecComp.USt(Itrans==3),pooledmap.UVStvecComp.VSt(Itrans==3),pooledmap.UVStvecComp.VecUSt(Itrans==3),pooledmap.UVStvecComp.VecVSt(Itrans==3),scale,'r');
    quiver(pooledmap.UVStvecComp.USt(Itrans==4),pooledmap.UVStvecComp.VSt(Itrans==4),pooledmap.UVStvecComp.VecUSt(Itrans==4),pooledmap.UVStvecComp.VecVSt(Itrans==4),scale,'k');
elseif COLORstate==1
    [Ms,Is]=min(angleUVSt,[],2);      % finds the canal that fits the data best
    quiver(pooledmap.UVStvecComp.USt(Is==1),pooledmap.UVStvecComp.VSt(Is==1),pooledmap.UVStvecComp.VecUSt(Is==1),pooledmap.UVStvecComp.VecVSt(Is==1),scale,'b');
    quiver(pooledmap.UVStvecComp.USt(Is==2),pooledmap.UVStvecComp.VSt(Is==2),pooledmap.UVStvecComp.VecUSt(Is==2),pooledmap.UVStvecComp.VecVSt(Is==2),scale,'m');
    quiver(pooledmap.UVStvecComp.USt(Is==3),pooledmap.UVStvecComp.VSt(Is==3),pooledmap.UVStvecComp.VecUSt(Is==3),pooledmap.UVStvecComp.VecVSt(Is==3),scale,'r');
elseif COLORstate==2
    angleUVStAll=[rad2deg(pooledmap.UVStvecComp.angleUVStASC),angleUVStInv(:,1),rad2deg(pooledmap.UVStvecComp.angleUVStPSC),...
        angleUVStInv(:,2),rad2deg(pooledmap.UVStvecComp.angleUVStLSC),angleUVStInv(:,3)];
    [Mall,Iall]=min(angleUVStAll,[],2);      % finds the canal that fits either the originl data or the inverse data best
    quiver(pooledmap.UVStvecComp.USt(Iall==1),pooledmap.UVStvecComp.VSt(Iall==1),pooledmap.UVStvecComp.VecUSt(Iall==1),pooledmap.UVStvecComp.VecVSt(Iall==1),scale,'b');
    quiver(pooledmap.UVStvecComp.USt(Iall==2),pooledmap.UVStvecComp.VSt(Iall==2),pooledmap.UVStvecComp.VecUSt(Iall==2),pooledmap.UVStvecComp.VecVSt(Iall==2),scale,'b');
    quiver(pooledmap.UVStvecComp.USt(Iall==3),pooledmap.UVStvecComp.VSt(Iall==3),pooledmap.UVStvecComp.VecUSt(Iall==3),pooledmap.UVStvecComp.VecVSt(Iall==3),scale,'m');
    quiver(pooledmap.UVStvecComp.USt(Iall==4),pooledmap.UVStvecComp.VSt(Iall==4),pooledmap.UVStvecComp.VecUSt(Iall==4),pooledmap.UVStvecComp.VecVSt(Iall==4),scale,'m');
    quiver(pooledmap.UVStvecComp.USt(Iall==5),pooledmap.UVStvecComp.VSt(Iall==5),pooledmap.UVStvecComp.VecUSt(Iall==5),pooledmap.UVStvecComp.VecVSt(Iall==5),scale,'r');
    quiver(pooledmap.UVStvecComp.USt(Iall==6),pooledmap.UVStvecComp.VSt(Iall==6),pooledmap.UVStvecComp.VecUSt(Iall==6),pooledmap.UVStvecComp.VecVSt(Iall==6),scale,'r');
elseif COLORstate==3
    quiver(pooledmap.UVStvecComp.USt,pooledmap.UVStvecComp.VSt,pooledmap.UVStvecComp.VecUSt,pooledmap.UVStvecComp.VecVSt,scale,'k','LineWidth',0.5);
end
set(f1,'position',[450 150 950 800]);
set(gcf, 'color', [1 1 1]);
%axis equal
xlim([-1 1]);ylim([-1 1]);
shg;
set(findobj(gcf, 'type','axes'), 'Visible','off')

try pooledmap.Itrans=Itrans; end    % the channel number (out of the four translation channels) that fits the data best

clusteringFlag=0;
if CLUSTERstate==2 && COLORstate==3
pooledmap=addClusterGridColor(pooledmap,COLORstate,f1,f2,f3,f4);
elseif CLUSTERstate==2 && COLORstate==0
pooledmap=addClusterGridColor(pooledmap,COLORstate,f1,f2,f3,f4);

elseif CLUSTERstate==1
pause(10);
clusteringFlag=get(handles.doneClusteringButton,'value');
while clusteringFlag==0
    pooledmap=addCluster(pooledmap,f1,f2,f3);
    pause(7);
    clusteringFlag=get(handles.doneClusteringButton,'value');
end

try
    % calculate the number of cells in each cluster that fit each channel
for z=1:size(pooledmap.clusterdata,2)       % number of clusters
    for t=1:4       % number of translation channels
cellCountPerChannel(z,t)=sum(pooledmap.clusterdata{1,z}.Itrans==t);
    end
end
sumCells=sum(cellCountPerChannel,2);

for g=1:size(cellCountPerChannel,1)
cellProportionPerChannel(g,:)=cellCountPerChannel(g,:)./sumCells(g,1);
end

pooledmap.cellCountPerChannel=cellCountPerChannel;
pooledmap.cellProportionPerChannel=cellProportionPerChannel;

figure;
f4=bar(pooledmap.cellProportionPerChannel)
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
        savefig(f1,['pooledMapAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_',transFieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        try saveas(f2,['pooledMapAlphaCorrClusters_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_',transFieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime],'fig'); end
        saveas(f3,['pooledMapAlphaCorrPolarHist_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_',transFieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        save(['pooledMapAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.mat'],'pooledmap');
    case 0
        saveas(f1,['pooledMap_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_',transFieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        try saveas(f2,['pooledMapClusters_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_',transFieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']); end
        saveas(f3,['pooledMapPolarHist_',num2str(filename),'_',num2str(discPoints),'_',fieldType,'_',transFieldType,'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.fig']);
        save(['pooledMap_',num2str(filename),'_',num2str(discPoints),'_DSI',num2str(DSIlow, '%0.2f'),'_',num2str(DSIhigh, '%0.2f'),'_created_',currenttime,'.mat'],'pooledmap');
end

%set(handles.alphaCorrbutton,'value',0);
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


% --- Executes on button press in IFbutton.
function IFbutton_Callback(hObject, eventdata, handles)
% hObject    handle to IFbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IFbutton


if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.IFbutton,'value',1);
else
     set(handles.IFbutton,'value',0);
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
     case 'iASCbutton'
         fieldType='iASC';
     case 'iPSCbutton'
         fieldType='iPSC';
     case 'iLSCbutton'
         fieldType='iLSC';
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
     case 'fourTransChannelColorbutton'
         COLORstate=0;
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
global DISTstate
global locDist
global widthDist
global mapType
global ABstate
global GSstate

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file or a combinedMap or a single-channel map file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file or a combinedMap file');
end

switch DATAMODELstate
    case 1
        histType=5;
    case 0
        histType=5;
    case 2
        histType=6;
    case 3
        histType=7;
end

calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,mapType,DISTstate,locDist,widthDist,ABstate,GSstate);


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
     case 'combinedMapbutton'
         DATAMODELstate=2;
         CLUSTNstate=0;
     case 'singleMapbutton'
         DATAMODELstate=3;
         CLUSTNstate=0;
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


% --- Executes on button press in fitSingleChannelMapsbutton.
function fitSingleChannelMapsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to fitSingleChannelMapsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2
global mapType
        
fitSurf2(polarity2,resolution,mapType);
        

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
global GSstate

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end

fieldType=1;
histType=4;
calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,fieldType,GSstate);


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
         set(handles.clusterNumber, 'enable', 'off');
     case 'bothPolbutton'
         polarity2=2;
         set(handles.clusterNumber, 'enable', 'on');
     case 'ASCLSCbutton'
         polarity2=3;
         set(handles.clusterNumber, 'enable', 'on');
     case 'realSingRotbutton'
         polarity2=4;
         set(handles.clusterNumber, 'enable', 'on');        
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
global TYPECLUSTstate
global TYPEFEATURECLUSTstate
global TYPESUPERCLUSTstate
global ONorONOFForOFFstate
global DSIlow
global DSIhigh
global NORMALRESstate
global postPcut
global magCut
global alphashift

% determine the current ONorONOFForOFFstate to be used for the clustering
% only.
ONorONOFForOFFstate=[];
if get(handles.ONDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,1]; end      % ONDS
if get(handles.ONOFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,2]; end   %ONOFFDS
if get(handles.OFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,3]; end     % OFFDS
if get(handles.ONtDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,4]; end     % ONtDS
if get(handles.ONsusOFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,5]; end     % ONsusOFFDS

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
        namesASC=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI0.17_1.00','_ASC'],0,'anywhere');    % DSI based on gratings
        for k=1:size(namesASC,2)
            ind=strfind(namesASC{k},'_');
            expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
            
            namesMap=getFileList(filedir,['vecCompAlphaCorr_',num2str(discPoints),'_DSI0.17_1.00'],0,'anywhere');      % DSI based on gratings 
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
            
            namesMap=getFileList(filedir,['vecComp_',num2str(discPoints),'_DSI0.17_1.00'],0,'anywhere');    % DSI based on gratings   
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
    
    
    % determine the minimum number of elements
elements=[];

for nu=1:numel(vecCompASC{k}.pResponse)
    elements=[elements,numel(vecCompASC{k}.pResponse{nu})];
end

% standardize the length of responses
for nu=1:numel(vecCompASC{k}.pResponse)
    numPresponseT{nu}=vecCompASC{k}.pResponse{nu}(1:min(elements));
end
DSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
clear numPresponseT

 maxRes=max(DSpResponse); 

 % a temporary fix for Hb9 and treatments
 DSIlow=0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    type=(vecCompASC{k}.DSI>=DSIlow) .* (vecCompASC{k}.DSI<=DSIhigh) .* (maxRes>=magCut);      % include all DS cells with DSI between DSI low and DSI high in the analysis
    
    for j=1:size(vecCompASC{k}.isCART,2)
        if logical(vecCompASC{k}.ONorONOFForOFF(j)) && logical(type(j))
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
        allCells.onset1(i,1)=vecCompASC{k}.onset1(j);
        
        end  
    end
end   % for analyzing all the retinas


%% unsupervised clustering on PCA coefficients
if TYPECLUSTstate  
% extract responses from all retinas

i=0;
%type=(allCells.DSI>=DSIlow) .* (allCells.DSI<=DSIhigh);      % include all DS cells with DSI between DSI low and DSI high in the analysis
    for j=1:size(allCells.DSI,1)    
            %if (sum(logical(allCells.ONorONOFForOFF(j)==ONorONOFForOFFstate)) && logical(type(j)))
        if (sum(logical(allCells.ONorONOFForOFF(j)==ONorONOFForOFFstate)))
            i=i+1;
            pResponse(i,1)=allCells.pResponse(j);
            ONorONOFForOFF(i,1)=allCells.ONorONOFForOFF(j);
            indpRes(i,1)=j;
        end
    end

%ONorONOFForOFF(size(ONorONOFForOFF,1)+1:(size(ONorONOFForOFF,1)+size(allCells.ONorONOFForOFF,1)-size(ONorONOFForOFF,1)),1)=0;    
    
plotAverageResponsePerType(ONorONOFForOFF,pResponse);
[idxT,inddataT]=clusteringDS4(pResponse,NORMALRESstate,postPcut);  
  
if sum(idxT==1)>sum(idxT==2)
    idxF(idxT==1)=2;
    idxF(idxT==2)=1;
else
    idxF=idxT;
end

if ~iscolumn(idxF); idxF=idxF'; end
idxF1(inddataT,1)=idxF;       % account for the cells that were omitted when applying the treshold for poterior probability. So, now, idxF1 corresponds to the same cells as in idxF.

% match is 1 if idxF2 and ONorONOFForOFF match, else it is zero.
ONorONOFForOFFt=ONorONOFForOFF;
idxF2=idxF1;
ONorONOFForOFFt(idxF2==0)=[];
idxF2(idxF2==0)=[];

match(idxF2==ONorONOFForOFFt)=1;
match(idxF2~=ONorONOFForOFFt)=0;
misclass=numel(match)-sum(match)       % number of misclassified cells
permisclass=100*(misclass./numel(match))   % percentage of misclassified cells

allCells.typeClust(indpRes,1)=idxF1;       % account for the missing cells in idxF. So, now, allCells.clustType corresponds to the same cells as in ONorONOFForOFF

end    

%% unsupervised clustering on extracted features
if TYPEFEATURECLUSTstate  
% extract responses from all retinas

i=0;
%type=(allCells.DSI>=DSIlow) .* (allCells.DSI<=DSIhigh);      % include all DS cells with DSI between DSI low and DSI high in the analysis
    for j=1:size(allCells.DSI,1)    
            %if (sum(logical(allCells.ONorONOFForOFF(j)==ONorONOFForOFFstate)) && logical(type(j)))
        if (sum(logical(allCells.ONorONOFForOFF(j)==ONorONOFForOFFstate)))
            i=i+1;
            pResponse(i,1)=allCells.pResponse(j);
            onset1(i,1)=allCells.onset1(j);
            ONorONOFForOFF(i,1)=allCells.ONorONOFForOFF(j);
            indpRes(i,1)=j;
        end
    end

%ONorONOFForOFF(size(ONorONOFForOFF,1)+1:(size(ONorONOFForOFF,1)+size(allCells.ONorONOFForOFF,1)-size(ONorONOFForOFF,1)),1)=0;    
    
plotAverageResponsePerType(ONorONOFForOFF,pResponse);

%[idxT,inddataT]=clusteringDS5(pResponse,onset1,NORMALRESstate,postPcut);  

[idxT,inddataT]=clusteringDS6(pResponse,onset1,NORMALRESstate,postPcut);  % for L-AP-4 experiment only
  
if sum(idxT==1)>sum(idxT==2)
    idxF(idxT==1)=2;
    idxF(idxT==2)=1;
else
    idxF=idxT;
end

if ~iscolumn(idxF); idxF=idxF'; end
idxF1(inddataT,1)=idxF;       % account for the cells that were omitted when applying the treshold for poterior probability. So, now, idxF1 corresponds to the same cells as in idxF.

% match is 1 if idxF2 and ONorONOFForOFF match, else it is zero.
ONorONOFForOFFt=ONorONOFForOFF;
idxF2=idxF1;
ONorONOFForOFFt(idxF2==0)=[];
idxF2(idxF2==0)=[];

match(idxF2==ONorONOFForOFFt)=1;
match(idxF2~=ONorONOFForOFFt)=0;
misclass=numel(match)-sum(match)       % number of misclassified cells
permisclass=100*(misclass./numel(match))   % percentage of misclassified cells

allCells.featureTypeClust(indpRes,1)=idxF1;       % account for the missing cells in idxF. So, now, allCells.clustType corresponds to the same cells as in ONorONOFForOFF

end    



%% supervised clustering
if TYPESUPERCLUSTstate  
% extract responses from all retinas

i=0;
%type=(allCells.DSI>=DSIlow) .* (allCells.DSI<=DSIhigh);      % include all DS cells with DSI between DSI low and DSI high in the analysis
    for j=1:size(allCells.DSI,1)    
        if (sum(logical(allCells.ONorONOFForOFF(j)==ONorONOFForOFFstate)))
            %if (sum(logical(allCells.ONorONOFForOFF(j)==ONorONOFForOFFstate)) && logical(type(j)))
            i=i+1;
            pResponse(i,1)=allCells.pResponse(j);
            ONorONOFForOFF(i,1)=allCells.ONorONOFForOFF(j);
            indpRes(i,1)=j;
        end
    end

%ONorONOFForOFF(size(ONorONOFForOFF,1)+1:(size(ONorONOFForOFF,1)+size(allCells.ONorONOFForOFF,1)-size(ONorONOFForOFF,1)),1)=0;    
    
plotAverageResponsePerType(ONorONOFForOFF,pResponse);


[superidxT,superinddataT]=superClusteringDS(pResponse,NORMALRESstate,postPcut);  
  
if sum(superidxT==1)>sum(superidxT==2)
    superidxF(superidxT==1)=2;
    superidxF(superidxT==2)=1;
else
    superidxF=superidxT;
end

if ~iscolumn(superidxF); superidxF=superidxF'; end
superidxF1(superinddataT,1)=superidxF;       % account for the cells that were omitted when applying the treshold for poterior probability. So, now, superidxF1 corresponds to the same cells as in superidxF.

% match is 1 if superidxF2 and ONorONOFForOFF match, else it is zero.
ONorONOFForOFFt=ONorONOFForOFF;
superidxF2=superidxF1;
ONorONOFForOFFt(superidxF2==0)=[];
superidxF2(superidxF2==0)=[];

match(superidxF2==ONorONOFForOFFt)=1;
match(superidxF2~=ONorONOFForOFFt)=0;
supermisclass=numel(match)-sum(match)       % number of misclassified cells
superpermisclass=100*(supermisclass./numel(match))   % percentage of misclassified cells

allCells.superTypeClust(indpRes,1)=superidxF1;       % account for the missing cells in superidxF. So, now, allCells.clustType corresponds to the same cells as in ONorONOFForOFF

end    

%%


filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');

switch ALPHACORRstate
    case 1
          save(['allCellsAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_PPC_',num2str(postPcut),'_DSI0.17_1.00','_created_',currenttime,'.mat'],'allCells');
    case 0
          save(['allCells_',num2str(filename),'_',num2str(discPoints),'_PPC_',num2str(postPcut),'_DSI0.17_1.00','_created_',currenttime,'.mat'],'allCells');
end


% --- Executes on button press in typeClustbutton.
function typeClustbutton_Callback(hObject, eventdata, handles)
% hObject    handle to typeClustbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of typeClustbutton

global TYPECLUSTstate

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.typeClustbutton,'value',1);
     TYPECLUSTstate=1;
else
     set(handles.typeClustbutton,'value',0);
     TYPECLUSTstate=0;
end


% --- Executes on button press in useTypeClustbutton.
function useTypeClustbutton_Callback(hObject, eventdata, handles)
% hObject    handle to useTypeClustbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useTypeClustbutton

global USETYPECLUSTstate

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.useTypeClustbutton,'value',1);
     USETYPECLUSTstate=1;
else
     set(handles.useTypeClustbutton,'value',0);
     USETYPECLUSTstate=0;
end


% --- Executes on button press in normalResbutton.
function normalResbutton_Callback(hObject, eventdata, handles)
% hObject    handle to normalResbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalResbutton

global NORMALRESstate

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.normalResbutton,'value',1);
     NORMALRESstate=1;
else
     set(handles.normalResbutton,'value',0);
     NORMALRESstate=0;
end


% --- Executes on button press in nullbutton.
function nullbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nullbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2
global numIter
global mapType

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end

histType=10;
calculateAllPossibleVecFieldsForNull(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2,mapType);
    
        



function postPcut_Callback(hObject, eventdata, handles)
% hObject    handle to postPcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of postPcut as text
%        str2double(get(hObject,'String')) returns contents of postPcut as a double

global postPcut

temp_postPcut=str2num(get(handles.postPcut,'String'));
if temp_postPcut<0 || temp_postPcut>1
    h=msgbox('cutoff must be between 0 and 1', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    postPcut=temp_postPcut;
end


% --- Executes during object creation, after setting all properties.
function postPcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to postPcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numIter_Callback(hObject, eventdata, handles)
% hObject    handle to numIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numIter as text
%        str2double(get(hObject,'String')) returns contents of numIter as a double
global numIter

temp_numIter=str2num(get(handles.numIter,'String'));
if temp_numIter<0 || temp_numIter>100000
    h=msgbox('number of iterations must be between 1 and 100000', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    numIter=temp_numIter;
end


% --- Executes during object creation, after setting all properties.
function numIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in singleChannelMapsbutton.
function singleChannelMapsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to singleChannelMapsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2
global mapType
global ABstate
global DISTstate
global locDist
global widthDist
global addNoise
global saculaPolarity 
global GSstate

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end

switch mapType
    case 1              % rotation

        % determine the source of alpha beta
        if ABstate==1          % get alpha beta from file
        [allPossibleFieldsName,PathName] = uigetfile('*.mat','Select a allPossibleFields rotation file of real data');
        load(allPossibleFieldsName);
        
        if GSstate==0
            alpharot=allPossibleFields.alphasing;
            betarot=allPossibleFields.betasing;
        elseif GSstate==1
            for g=1:size(allPossibleFields.alphasing,2)
                
                if  allPossibleFields.betasing(g)==0
                    allPossibleFields.betasing(g)=180;
                else
                    allPossibleFields.betasing(g)=mod(180-allPossibleFields.betasing(g),180);
                end
                
                [alphaRG,betaRG]=rotateNormalToExtraPersonalSpace3('Right',allPossibleFields.alphasing(g),allPossibleFields.betasing(g));
                [alphaLG,betaLG]=rotateNormalToExtraPersonalSpace3('Left',allPossibleFields.alphasing(g),allPossibleFields.betasing(g));
                alpharot(1,g)=alphaRG;
                betarot(1,g)=betaRG;
                alpharotL(1,g)=alphaLG;
                betarotL(1,g)=betaLG;
            end
        end
        
        
        histType=20;
        fieldName=[];
        
        alpha=alpharot(1);
        beta=betarot(1);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)

        alpha=alpharot(2);
        beta=betarot(2);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)

        alpha=alpharot(3);
        beta=betarot(3);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)

        alpha=alpharot(4);
        beta=betarot(4);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        try
        alpha=alpharot(5);
        beta=betarot(5);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)

        alpha=alpharot(6);
        beta=betarot(6);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
         
        alpha=alpharotL(1);
        beta=betarotL(1);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=alpharotL(2);
        beta=betarotL(2);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=alpharotL(3);
        beta=betarotL(3);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=alpharotL(4);
        beta=betarotL(4);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=alpharotL(5);
        beta=betarotL(5);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=alpharotL(6);
        beta=betarotL(6);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        end
        
        elseif ABstate==2          % use alpha beta of real canals
              
switch polarity2
    case 1
        histType=20;
        fieldName='ASC';
        %normal=[0.596,0.766,0.241];         %ASC
        normal=[0.604,0.758,0.243]; % ASC symmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,fieldType,fieldName,alpha,beta,GSstate)
        
        fieldName='PSC';
        %normal=[0.471,-0.766,0.438];        %PSC
        normal=[0.447,-0.76,0.469];    % PSC symmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,fieldType,fieldName,alpha,beta,GSstate)
        
        fieldName='LSC';
        %normal=[0.421,-0.060,-0.905];       %LSC
        normal=[0.449,-0.085,-0.886];   % LSC symmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,fieldType,fieldName,alpha,beta,GSstate)
        
    case 2
        
        histType=20;
        fieldName='ASC';
        
        %normal=[0.596,0.766,0.241];         %ASC
        normal=[0.604,0.758,0.243]; % ASC ymmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
                                     
        fieldName='PSC';
        %normal=[0.471,-0.766,0.438];        %PSC
        normal=[0.447,-0.76,0.469];    % PSC symmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        fieldName='LSC';
        %normal=[0.421,-0.060,-0.905];       %LSC
        normal=[0.449,-0.085,-0.886];   % LSC symmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        fieldName='iASC';
        %normal=[0.596,0.766,0.241];         %iASC
        normal=[0.604,0.758,0.243]; % iASC ymmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        beta=90-abs(90-beta);
        if alpha<=180
            alpha=180+alpha;
        else alpha>180 
            alpha=alpha-180;
        end
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        fieldName='iPSC';
        %normal=[0.471,-0.766,0.438];          %iPSC
        normal=[0.447,-0.76,0.469];    % iPSC symmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        beta=90-abs(90-beta);
        if alpha<=180
            alpha=180+alpha;
        else alpha>180 
            alpha=alpha-180;
        end
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        fieldName='iLSC';
        %normal=[0.421,-0.060,-0.905];          %iLSC
        normal=[0.449,-0.085,-0.886];   % iLSC symmetric
        [alpharot,betarot]=AB(normal);
        alpha=rad2deg(alpharot);
        beta=rad2deg(betarot);
        beta=90-abs(90-beta);
        if alpha<=180
            alpha=180+alpha;
        else alpha>180 
            alpha=alpha-180;
        end
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)

end
        end

    case 0              % translation
        
        % determine the source of alpha beta
        if ABstate==1          % get alpha beta from allPossibleFields file
            [allPossibleFieldsName,PathName] = uigetfile('*.mat','Select a allPossibleFields translation file of real data');
            load(allPossibleFieldsName);
            
            if GSstate==0
                alphatrans=allPossibleFields.alphasing;
                betatrans=allPossibleFields.betasing;
            elseif GSstate==1
                for g=1:size(allPossibleFields.alphasing,2)
                    
                if  allPossibleFields.betasing(g)==0
                    allPossibleFields.betasing(g)=180;
                else
                    allPossibleFields.betasing(g)=mod(180-allPossibleFields.betasing(g),180);
                end
                    
                    [alphaRG,betaRG]=rotateNormalToExtraPersonalSpace3('Right',allPossibleFields.alphasing(g),allPossibleFields.betasing(g));
                    [alphaLG,betaLG]=rotateNormalToExtraPersonalSpace3('Left',allPossibleFields.alphasing(g),allPossibleFields.betasing(g));
                    alphatrans(1,g)=alphaRG;
                    betatrans(1,g)=betaRG;
                    alphatransL(1,g)=alphaLG;
                    betatransL(1,g)=betaLG;
                end
            end
 
        elseif ABstate==0          % use fixed alpha beta
            alphatrans=[0,90,180,270];
            betatrans=[90,90,90,90];    
        elseif ABstate==3          % use alpha beta of sacula     
            switch saculaPolarity
                case 1      
                    %normal=[-0.144,-0.315,0.938];       % vector of the sacular macula
                    normal=[-0.119,-0.286,0.951];       % vector of the sacular macula Symmetric
                    [alpharot,betarot]=AB(normal);
                    alphaT=rad2deg(alpharot);
                    beta=rad2deg(betarot);
                case -1
                    %normal=-[-0.144,-0.315,0.938];       % vector of the inverse sacular macula
                    normal=-[-0.119,-0.286,0.951];       % vector of the inverse sacular macula Symmetric
                    [alpharot,betarot]=AB(normal);
                    alphaT=rad2deg(alpharot);
                    beta=rad2deg(betarot);
                    beta=90-abs(90-beta);
                    if alphaT<=180
                        alphaT=180+alphaT;
                    else alphaT>180
                        alphaT=alphaT-180;
                    end
            end
        
            alphatrans=floor(alphaT);
            betatrans=floor(beta);
            
        end
        
             
        histType=20;
        fieldName=[];
        
        alpha=floor(alphatrans(1));
        beta=floor(betatrans(1));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        try
        alpha=floor(alphatrans(2));
        beta=floor(betatrans(2));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)

        alpha=floor(alphatrans(3));
        beta=floor(betatrans(3));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)

        alpha=floor(alphatrans(4));
        beta=floor(betatrans(4));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=floor(alphatrans(5));
        beta=floor(betatrans(5));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=floor(alphatrans(6));
        beta=floor(betatrans(6));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=floor(alphatransL(1));
        beta=floor(betatransL(1));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=floor(alphatransL(2));
        beta=floor(betatransL(2));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=floor(alphatransL(3));
        beta=floor(betatransL(3));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        
        alpha=floor(alphatransL(4));
        beta=floor(betatransL(4));
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alpha,beta,GSstate)
        end
end

% --- Executes on button press in stepwisebutton.
function stepwisebutton_Callback(hObject, eventdata, handles)
% hObject    handle to stepwisebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2
global numIter
global addNoise

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end

switch polarity2
    case 1
        histType=20;
        fieldType=1;        % rotational vector field
        calculateVecFieldStepwiseWithNoise(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2,addNoise,fieldType);
    case 2
        
        histType=20;
        fieldType=1;        % rotational vector field
        calculateVecFieldStepwiseWithNoise(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2,addNoise,fieldType);
end



function addNoise_Callback(hObject, eventdata, handles)
% hObject    handle to addNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of addNoise as text
%        str2double(get(hObject,'String')) returns contents of addNoise as a double
global addNoise

temp_addNoise=str2num(get(handles.addNoise,'String'));
if temp_addNoise<0 || temp_addNoise>360
    h=msgbox('standard deviation of noise must be between 0 and 360', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    addNoise=temp_addNoise;
end


% --- Executes during object creation, after setting all properties.
function addNoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to addNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stepwiseTransbutton.
function stepwiseTransbutton_Callback(hObject, eventdata, handles)
% hObject    handle to stepwiseTransbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2
global numIter
global addNoise

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end

switch polarity2
    case 1
        histType=20;
        fieldType=0;        % translational vector field
        calculateVecFieldStepwiseWithNoise(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2,addNoise,fieldType);
    case 2
        
        histType=20;
        fieldType=0;        % translational vector field
        calculateVecFieldStepwiseWithNoise(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2,addNoise,fieldType);
end




% --- Executes on button press in trans4Channelbutton.
function trans4Channelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to trans4Channelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2
global GSstate

switch ALPHACORRstate
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

    case 0
        [FileName,PathName] = uigetfile('*.mat','Select a pooledMap file');
end


        histType=20;
        fieldType=0;
        alpha=0;
        beta=90;
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,fieldType,alpha,beta,GSstate)
        
        histType=20;
        fieldType=0;
        alpha=90;
        beta=90;
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,fieldType,alpha,beta,GSstate)
        
        histType=20;
        fieldType=0;
        alpha=180;
        beta=90;
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,fieldType,alpha,beta,GSstate)
        
        histType=20;
        fieldType=0;
        alpha=270;
        beta=90;
        calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,fieldType,alpha,beta,GSstate)
        
        


% --- Executes on button press in typeSuperClustbutton.
function typeSuperClustbutton_Callback(hObject, eventdata, handles)
% hObject    handle to typeSuperClustbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of typeSuperClustbutton

global TYPESUPERCLUSTstate

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.typeSuperClustbutton,'value',1);
     TYPESUPERCLUSTstate=1;
else
     set(handles.typeSuperClustbutton,'value',0);
     TYPESUPERCLUSTstate=0;
end


% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global USETYPECLUSTstate

SUPERVISIONButton=get(eventdata.NewValue,'tag');
switch SUPERVISIONButton
     case 'nonSupervisedbutton'
         USETYPECLUSTstate=1;
     case 'supervisedbutton'
         USETYPECLUSTstate=2;
     case 'manualbutton'
         USETYPECLUSTstate=0;
     case 'featureNonSupervisedbutton'
         USETYPECLUSTstate=3;     
end



function magCut_Callback(hObject, eventdata, handles)
% hObject    handle to magCut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magCut as text
%        str2double(get(hObject,'String')) returns contents of magCut as a double

global magCut

temp_magCut=str2num(get(handles.magCut,'String'));
if temp_magCut<0
    h=msgbox('magnitude threshold must be larger than 0', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    magCut=temp_magCut;
end


% --- Executes during object creation, after setting all properties.
function magCut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magCut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in typeFeatureClustbutton.
function typeFeatureClustbutton_Callback(hObject, eventdata, handles)
% hObject    handle to typeFeatureClustbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of typeFeatureClustbutton

global TYPEFEATURECLUSTstate

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.typeFeatureClustbutton,'value',1);
     TYPEFEATURECLUSTstate=1;
else
     set(handles.typeFeatureClustbutton,'value',0);
     TYPEFEATURECLUSTstate=0;
end


% --- Executes when selected object is changed in uibuttongroup4.
function uibuttongroup4_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup4 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global transFieldType

transFieldTypeButton=get(eventdata.NewValue,'tag');
switch transFieldTypeButton
     case 'nasalbutton'
         transFieldType='nasal';
     case 'dorsalbutton'
         transFieldType='dorsal';
     case 'temporalbutton'
         transFieldType='temporal';
     case 'ventralbutton'
         transFieldType='ventral';
     case 'allTransbutton'
         transFieldType='all';
     case 'nonTransbutton'
         transFieldType=[];
end



function alphashift_Callback(hObject, eventdata, handles)
% hObject    handle to alphashift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphashift as text
%        str2double(get(hObject,'String')) returns contents of alphashift as a double
global alphashift

temp_alphashift=str2num(get(handles.alphashift,'String'));
if temp_alphashift<0 || temp_alphashift>360
    h=msgbox('shift must be between 0 and 360 degrees', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    alphashift=temp_alphashift;
end


% --- Executes during object creation, after setting all properties.
function alphashift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphashift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup5.
function uibuttongroup5_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup5 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ABstate

alphabetaSourceButton=get(eventdata.NewValue,'tag');
switch alphabetaSourceButton
     case 'alphabetaFixedbutton'
         ABstate=0;
         set(handles.alphashift, 'enable', 'on')
         set(handles.saculaPolarity, 'enable', 'off')
     case 'alphabetaFromFilebutton'
         ABstate=1;
         set(handles.alphashift, 'enable', 'off')
         set(handles.saculaPolarity, 'enable', 'off')
     case 'realCanalsbutton'
         ABstate=2;
         set(handles.alphashift, 'enable', 'off')
         set(handles.saculaPolarity, 'enable', 'off')
     case 'alphabetSaculabutton'    
         ABstate=3;
         set(handles.alphashift, 'enable', 'off')
         set(handles.saculaPolarity, 'enable', 'on')
end       



function locDist_Callback(hObject, eventdata, handles)
% hObject    handle to locDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of locDist as text
%        str2double(get(hObject,'String')) returns contents of locDist as a double
global locDist

temp_locDist=str2num(get(handles.locDist,'String'));
if temp_locDist<0 || temp_locDist>360
    h=msgbox('location of cutoff must be between 0 and 360 degrees', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    locDist=temp_locDist;
end


% --- Executes during object creation, after setting all properties.
function locDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to locDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function widthDist_Callback(hObject, eventdata, handles)
% hObject    handle to widthDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of widthDist as text
%        str2double(get(hObject,'String')) returns contents of widthDist as a double
global widthDist

temp_widthDist=str2num(get(handles.widthDist,'String'));
if temp_widthDist<0 || temp_widthDist>360
    h=msgbox('width of cutoff must be between 0 and 360 degrees', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    widthDist=temp_widthDist;
end


% --- Executes during object creation, after setting all properties.
function widthDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to widthDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup6.
function uibuttongroup6_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup6 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DISTstate

distanceButton=get(eventdata.NewValue,'tag');
switch distanceButton
     case 'noDistFuncbutton'
         DISTstate=0;
         set(handles.locDist, 'enable', 'off')
         set(handles.widthDist, 'enable', 'off')
     case 'applyDistFuncbutton'
         DISTstate=1;
         set(handles.locDist, 'enable', 'on')
         set(handles.widthDist, 'enable', 'on')
         
         
end    


% --- Executes on button press in combineSingleMapsbutton.
function combineSingleMapsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to combineSingleMapsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALPHACORRstate
global DATAMODELstate
global CLUSTNstate
global resolution
global polarity2
global mapType
global DISTstate
global locDist
global widthDist
global ABstate


combineSingleMaps(polarity2,mapType,DISTstate,locDist,widthDist,ABstate)


% --- Executes when selected object is changed in uibuttongroup7.
function uibuttongroup7_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup7 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mapType

mapTypeButton=get(eventdata.NewValue,'tag');
switch mapTypeButton
     case 'rotbutton'
         mapType=1;    
     case 'transbutton'
         mapType=0;  
     case 'mixedTransRotbutton'
         mapType=2;
end    



function saculaPolarity_Callback(hObject, eventdata, handles)
% hObject    handle to saculaPolarity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saculaPolarity as text
%        str2double(get(hObject,'String')) returns contents of saculaPolarity as a double
global saculaPolarity 

temp_saculaPolarity=str2num(get(handles.saculaPolarity,'String'));
if temp_saculaPolarity<-1 || temp_saculaPolarity>1
    h=msgbox('polarity must be either 1 or -1', 'Error','error');
else
    %set(handles.Alpha,'String', num2str(temp_alpha));
    saculaPolarity=temp_saculaPolarity;
end

% --- Executes during object creation, after setting all properties.
function saculaPolarity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saculaPolarity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in globalSpacebutton.
function globalSpacebutton_Callback(hObject, eventdata, handles)
% hObject    handle to globalSpacebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of globalSpacebutton
global GSstate

if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.globalSpacebutton,'value',1);
     GSstate=1;
else
     set(handles.globalSpacebutton,'value',0);
     GSstate=0;
end
