function varargout = clusteringGUI(varargin)
% CLUSTERINGGUI MATLAB code for clusteringGUI.fig
%      CLUSTERINGGUI, by itself, creates a new CLUSTERINGGUI or raises the existing
%      singleton*.
%
%      H = CLUSTERINGGUI returns the handle to a new CLUSTERINGGUI or the handle to
%      the existing singleton*.
%
%      CLUSTERINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTERINGGUI.M with the given input arguments.
%
%      CLUSTERINGGUI('Property','Value',...) creates a new CLUSTERINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before clusteringGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to clusteringGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help clusteringGUI

% Last Modified by GUIDE v2.5 07-Mar-2015 14:23:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @clusteringGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @clusteringGUI_OutputFcn, ...
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


% --- Executes just before clusteringGUI is made visible.
function clusteringGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to clusteringGUI (see VARARGIN)

% Choose default command line output for clusteringGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes clusteringGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%% initialize_input_parameters
clear global"
global cellTypestate
global ONorONOFForOFFstate

%initialize handles and display default parameters on user interface
set(handles.DSbutton,'value',1);
cellTypestate=1;
set(handles.ONDSbutton,'value',1);
ONorONOFForOFFstate=1;


% --- Outputs from this function are returned to the command line.
function varargout = clusteringGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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




% --- Executes on button press in StartCluteringbutton.
function StartCluteringbutton_Callback(hObject, eventdata, handles)
% hObject    handle to StartCluteringbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global cellTypestate
global ONorONOFForOFFstate


% determine the current Cell type state
cellTypestate=[];
if get(handles.DSbutton,'value')==1 cellTypestate=[cellTypestate,1]; end      % DS
if get(handles.OSbutton,'value')==1 cellTypestate=[cellTypestate,2]; end   %OS
if get(handles.nonDSnonOSbutton,'value')==1 cellTypestate=[cellTypestate,3]; end     % nonDSnonOS

% determine the current ONorONOFForOFFstate
ONorONOFForOFFstate=[];
if get(handles.ONDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,1]; end      % ONDS
if get(handles.ONOFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,2]; end   %ONOFFDS
if get(handles.OFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,3]; end     % OFFDS
if get(handles.ONtDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,4]; end     % ONtDS
if get(handles.ONsusOFFDSbutton,'value')==1 ONorONOFForOFFstate=[ONorONOFForOFFstate,5]; end     % ONsusDS

filedir=cd;
namesMap=getFileList(filedir,'map_',0,'anywhere');
 
% extract responses from all retinas
x=0;
i=0;
for k=1:size(namesMap,2)       % number of retinas
    load(num2str(namesMap{k}));
    for j=1:size(map.isDS,1)
        if map.isDS(j)==1; cellType(j)=1; end
        if map.isOS(j)==1; cellType(j)=2; end
        if map.isDS(j)==0 && map.isOS(j)==0; cellType(j)=3; end
        if (sum(logical(map.ONorONOFForOFF(j)==ONorONOFForOFFstate)) && sum(logical(cellType(j)==cellTypestate)))
            i=i+1;
            pResponse(i,1)=map.Presponse(j);
            ONorONOFForOFF(i,1)=map.ONorONOFForOFF(j);
        elseif cellTypestate==[1,2,3]
            i=i+1;
            pResponse(i,1)=map.Presponse(j);
            ONorONOFForOFF(i,1)=map.ONorONOFForOFF(j);
        end
    end
end

plotAverageResponsePerType(ONorONOFForOFF,pResponse);
[idx]=clusteringDS4(pResponse);   


if sum(idx==1)>sum(idx==2)
    idxT(idx==1)=2;
    idxT(idx==2)=1;
else
    idxT=idx;
end

% match is 1 if idx and ONorONOFForOFF match, else it is zero.
match(idxT==ONorONOFForOFF)=1;
match(idxT~=ONorONOFForOFF)=0;
misclass=numel(match)-sum(match);       % number of misclassified cells
permisclass=100*(misclass./numel(match));   % percentage of misclassified cells



     






    
    
   



% --- Executes on button press in DSbutton.
function DSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to DSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DSbutton
if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.DSbutton,'value',1);
else
     set(handles.DSbutton,'value',0)
end


% --- Executes on button press in OSbutton.
function OSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OSbutton
if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.OSbutton,'value',1);
else
     set(handles.OSbutton,'value',0)
end


% --- Executes on button press in nonDSnonOSbutton.
function nonDSnonOSbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nonDSnonOSbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonDSnonOSbutton
if (get(hObject,'Value') == get(hObject,'Max'))
	 set(handles.nonDSnonOSbutton,'value',1);
else
     set(handles.nonDSnonOSbutton,'value',0)
end




