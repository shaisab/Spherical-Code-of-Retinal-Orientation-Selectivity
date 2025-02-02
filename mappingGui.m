function varargout = mappingGui(varargin)
% MAPPINGGUI MATLAB code for mappingGui.fig
%      MAPPINGGUI, by itself, creates a new MAPPINGGUI or raises the existing
%      singleton*.
%
%      H = MAPPINGGUI returns the handle to a new MAPPINGGUI or the handle to
%      the existing singleton*.
%
%      MAPPINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPPINGGUI.M with the given input arguments.
%
%      MAPPINGGUI('Property','Value',...) creates a new MAPPINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mappingGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mappingGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mappingGui

% Last Modified by GUIDE v2.5 12-Feb-2016 18:19:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mappingGui_OpeningFcn, ...
    'gui_OutputFcn',  @mappingGui_OutputFcn, ...
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


% --- Executes just before mappingGui is made visible.
function mappingGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mappingGui (see VARARGIN)

% Choose default command line output for mappingGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mappingGui wait for user response (see UIRESUME)
% uiwait(handles.mappingGui);

%%%%%% initialize_input_parameters
global alpha
global beta
global discPoints
global visAngle
global polarity
global rotateFlag
global tangentFlag
global radialFlag
global fieldType
global cellType
global DSIhigh
global DSIlow

discPoints=40;
visAngle=180;
DSIhigh=1;
DSIlow=0.17;
normal=[0.592,0.764,0.247];     % initialize alpha and beta for ASC (when no button is pressed)
%normal=[0.724,0.615,0.300];     % initialize alpha and beta for ASC Prime (when no button is pressed)
[alpha,beta]=AB(normal);
polarity=1;
rotateFlag=1;
tangentFlag=0;
radialFlag=0;
cellType='ONDS';
fieldType='ASC';

%display default parameters on user interface
set(handles.discPoints,'String', num2str(discPoints));

%display default parameters on user interface
set(handles.visAngle,'String', num2str(visAngle));

%display default parameters on user interface
set(handles.DSIhigh,'String', num2str(DSIhigh));

%display default parameters on user interface
set(handles.DSIlow,'String', num2str(DSIlow));





% --- Outputs from this function are returned to the command line.
function varargout = mappingGui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in retinalMapbutton.
function retinalMapbutton_Callback(hObject, eventdata, handles)
% hObject    handle to retinalMapbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mapCellsCorr_ss2;


% --- Executes on button press in retinalMapAlphaCorrbutton.
function retinalMapAlphaCorrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to retinalMapAlphaCorrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mapCellsAlphaCorr_HUJI;


% --- Executes on button press in mappingbutton.
function mappingbutton_Callback(hObject, eventdata, handles)
% hObject    handle to mappingbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%global discPoints

discPoints=str2num(get(handles.discPoints,'String'));

Mapping(discPoints);


% --- Executes on button press in mappingAlphaCorrbutton.
function mappingAlphaCorrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to mappingAlphaCorrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
discPoints=str2num(get(handles.discPoints,'String'));

MappingAlphaCorr(discPoints);


% --- Executes on button press in vecorFieldbutton.
function vecorFieldbutton_Callback(hObject, eventdata, handles)
% hObject    handle to vecorFieldbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global alpha
global beta
global discPoints
global visAngle
global polarity
global rotateFlag
global tangentFlag
global radialFlag
global fieldType

if radialFlag
    VecDataRad(polarity, discPoints);
else
    VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints, fieldType,visAngle);
end
rotateFlag=1;
tangentFlag=0;
radialFlag=0;





% --- Executes on button press in vecComparisonbutton.
function vecComparisonbutton_Callback(hObject, eventdata, handles)
% hObject    handle to vecComparisonbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cellType
global fieldType

discPoints=str2num(get(handles.discPoints,'String'));

% if ~exist('cellType','var')
%     cellType='ONDS';
% end

VecComparison(discPoints,cellType,fieldType)


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global cellType

cellTypeButton=get(eventdata.NewValue,'tag');
switch cellTypeButton
    case 'ONDSbutton'
        cellType='ONDS';
    case 'ONOFFDSbutton'
        cellType='ONOFFDS';
    case 'retroONDSbutton'
        cellType='retroONDS';
end



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
    %set(handles.Alpha,'String', num2str(temp_alpha));
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


% --- Executes on button press in AllcardinalONOFFbutton.
function AllcardinalONOFFbutton_Callback(hObject, eventdata, handles)
% hObject    handle to AllcardinalONOFFbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global alpha
global beta
global discPoints
global polarity
global rotateFlag
global tangentFlag
global radialFlag

error=0;
fixedONOFFdown=0+error;

alpha=fixedONOFFdown+pi;             %fixed ON-OFF up
beta=pi/2;
rotateFlag=1;
VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,'cardinalUp');
VecComparisonAll('ONOFFDS','cardinalUp');

alpha=fixedONOFFdown+(pi+pi/2);          %fixed ON-OFF left
beta=pi/2;
rotateFlag=1;
VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,'cardinalLeft');
VecComparisonAll('ONOFFDS','cardinalLeft');

alpha=fixedONOFFdown;                 %fixed ON-OFF down
beta=pi/2;
rotateFlag=1;
VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,'cardinalDown');
VecComparisonAll('ONOFFDS','cardinalDown');

alpha=fixedONOFFdown+(pi/2);       %fixed ON-OFF right
beta=pi/2;
rotateFlag=1;
VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,'cardinalRight');
VecComparisonAll('ONOFFDS','cardinalRight');

% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

global alpha
global beta
global polarity
global rotateFlag
global tangentFlag
global radialFlag
global fieldType


% setup parameters for fixed preferred direction models
error=0;
fixedLSC=-1.7+error;
fixedONOFFdown=0+error;

% selection of models
fieldTypeButton=get(eventdata.NewValue,'tag');
switch fieldTypeButton
    case 'ASCbutton'
        normal=[0.592,0.764,0.247];     %ASC
        %normal=[0.724,0.615,0.300];      %ASC Prime
        [alpha,beta]=AB(normal);
        rotateFlag=1;
        fieldType='ASC';
    case 'PSCbutton'
        normal=[0.470,-0.764,0.438];       %PSC
        %normal=[0.678,-0.638,0.359];        %PSC Prime
        [alpha,beta]=AB(normal);
        rotateFlag=1;
        fieldType='PSC';
    case 'LSCbutton'
        normal=[0.421,-0.056,-0.901];       %LSC
        %normal=[0.533,-0.102,-0.829];        %LSC Prime
        [alpha,beta]=AB(normal);
        rotateFlag=1;
        fieldType='LSC';
    case 'fixedupbutton'
        alpha=-deg2rad(120)+fixedLSC;         %fixed up-forward
        beta=pi/2;
        rotateFlag=1;
        fieldType='fixedUpForward';
    case 'fixeddownbutton'
        alpha=deg2rad(120)+fixedLSC;          %fixed down-backward
        beta=pi/2;
        rotateFlag=1;
        fieldType='fixedDownBackwards';
    case 'fixedbackwardsbutton'
        alpha=fixedLSC;                       %fixed horizontal
        beta=pi/2;
        rotateFlag=1;
        fieldType='fixedHorizontal';
    case 'tangentCWbutton'
        alpha=0;                              % tanget CW
        beta=0;
        polarity=-1;
        tangentFlag=1;
        fieldType='tangentCW';
    case 'tangentCCWbutton'
        alpha=0;                            % tangent CCW
        beta=0;
        polarity=1;
        tangentFlag=1;
        fieldType='tangentCCW';
    case 'radialINbutton'
        polarity=-1;                         % radial inward
        radialFlag=1;
        fieldType='radialInward';
    case 'radialOUTbutton'
        polarity=1;                          % radial outward
        radialFlag=1;
        fieldType='radialOutward';
    case 'upbutton'
        alpha=fixedONOFFdown+pi;             %fixed ON-OFF up
        beta=pi/2;
        rotateFlag=1;
        fieldType='fixedONOFFup';
    case 'leftbutton'
        alpha=fixedONOFFdown+(pi+pi/2);          %fixed ON-OFF left
        beta=pi/2;
        rotateFlag=1;
        fieldType='fixedONOFFleft';
    case 'downbutton'
        alpha=fixedONOFFdown;                 %fixed ON-OFF down
        beta=pi/2;
        rotateFlag=1;
        fieldType='fixedONOFFdown';
    case 'rightbutton'
        alpha=fixedONOFFdown+(pi/2);       %fixed ON-OFF right
        beta=pi/2;
        rotateFlag=1;
        fieldType='fixedONOFFright';
end


% --- Executes on button press in standardMappingbutton.
function standardMappingbutton_Callback(hObject, eventdata, handles)
% hObject    handle to standardMappingbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

discPoints=str2num(get(handles.discPoints,'String'));

standardMapping(discPoints);


% --- Executes on button press in canalsbutton.
function canalsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to canalsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global alpha
global beta
global discPoints
global visAngle
global polarity
global rotateFlag
global DSIhigh
global DSIlow

%error=deg2rad(45);
error=0;

normal=[0.592,0.764,0.247];     %ASC
%normal=[0.724,0.615,0.300];      %ASC Prime
[alpha,beta]=AB(normal);
alpha=alpha+error;
rotateFlag=1;
polarity=1;
VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,'ASC',visAngle);
VecComparisonInterp(discPoints,'ASC',DSIlow,DSIhigh);

normal=[0.470,-0.764,0.438];       %PSC
%normal=[0.678,-0.638,0.359];        %PSC Prime
[alpha,beta]=AB(normal);
alpha=alpha+error;
rotateFlag=1;
polarity=1;
VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,'PSC',visAngle);
VecComparisonInterp(discPoints,'PSC',DSIlow,DSIhigh);

normal=[0.421,-0.056,-0.901];       %LSC
%normal=[0.533,-0.102,-0.829];       %LSC Prime
[alpha,beta]=AB(normal);
alpha=alpha+error;
rotateFlag=1;
polarity=1;
VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,'LSC',visAngle);
VecComparisonInterp(discPoints,'LSC',DSIlow,DSIhigh);


filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
%filedate=filedir(ind(2)+1:ind(3)-1);

% canals={'ASC','PSC','LSC'};
% type={'ONDS','ONOFFDS'};
% %type={'ONDS','ONOFFDS','retro'};
% for j=1:size(type,2)
%     for i=1:size(canals,2)
%         names=getFileList(filedir,['vecComp_',num2str(discPoints),'_',type{j},'_',canals{i},'_'],0,'anywhere');
%         thename=shortfile(cell2mat(names));
%         thename=thename(1:strfind(thename,filename)-2);
%         vecComp.(thename)=load(cell2mat(names));
%         delete(cell2mat(names));
%     end
% end
%     thename=thename(1:strfind(thename,type{j})-2);
%     save([thename,'_',filename,'.mat'],'vecComp');
    
    close all;



    
    % --- Executes on button press in canalsAlphaCorrbutton.
function canalsAlphaCorrbutton_Callback(hObject, eventdata, handles)
% hObject    handle to canalsAlphaCorrbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global alpha
global beta
global discPoints
global visAngle
global polarity
global rotateFlag
global DSIhigh
global DSIlow

error=0;

normal=[0.592,0.764,0.247];     %ASC
%normal=[0.724,0.615,0.300];    %ASC Prime
[alpha,beta]=AB(normal);
alpha=alpha+error;
rotateFlag=1;
polarity=1;
VectorFieldGeneralizedForGuiAlphaCorr(alpha, beta, polarity, discPoints,'ASC',visAngle);
VecComparisonInterpAlphaCorr(discPoints,'ASC',DSIlow,DSIhigh);

normal=[0.470,-0.764,0.438];       %PSC
%normal=[0.678,-0.638,0.359];      %PSC Prime
[alpha,beta]=AB(normal);
alpha=alpha+error;
rotateFlag=1;
polarity=1;
VectorFieldGeneralizedForGuiAlphaCorr(alpha, beta, polarity, discPoints,'PSC',visAngle);
VecComparisonInterpAlphaCorr(discPoints,'PSC',DSIlow,DSIhigh);

normal=[0.421,-0.056,-0.901];       %LSC
%normal=[0.533,-0.102,-0.829];      %LSC Prime
[alpha,beta]=AB(normal);
alpha=alpha+error;
rotateFlag=1;
polarity=1;
VectorFieldGeneralizedForGuiAlphaCorr(alpha, beta, polarity, discPoints,'LSC',visAngle);
VecComparisonInterpAlphaCorr(discPoints,'LSC',DSIlow,DSIhigh);


filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
%filedate=filedir(ind(2)+1:ind(3)-1);
    close all;






    
    
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


% --- Executes on button press in alphaCorrConvolutionbutton.
function alphaCorrConvolutionbutton_Callback(hObject, eventdata, handles)
% hObject    handle to alphaCorrConvolutionbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

alphaCorr4;
