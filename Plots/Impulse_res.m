function varargout = Impulse_res(varargin)
% IMPULSE_RES MATLAB code for Impulse_res.fig
%      IMPULSE_RES, by itself, creates a new IMPULSE_RES or raises the existing
%      singleton*.
%
%      H = IMPULSE_RES returns the handle to a new IMPULSE_RES or the handle to
%      the existing singleton*.
%
%      IMPULSE_RES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPULSE_RES.M with the given input arguments.
%
%      IMPULSE_RES('Property','Value',...) creates a new IMPULSE_RES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Impulse_res_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Impulse_res_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Impulse_res

% Last Modified by GUIDE v2.5 19-Feb-2015 18:50:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Impulse_res_OpeningFcn, ...
                   'gui_OutputFcn',  @Impulse_res_OutputFcn, ...
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


% --- Executes just before Impulse_res is made visible.
function Impulse_res_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Impulse_res (see VARARGIN)

% Choose default command line output for Impulse_res
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Impulse_res wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global ii;
ii=1;
global h_T;
global t_sample;
global timei;
timei=0:t_sample:size(h_T,1)*t_sample;
global d_t;
global tauv;
tauv=0:d_t*1e6:(size(h_T,2)-1)*d_t*1e6;

set(handles.t_sample,'String',t_sample*1e3);
set(handles.fl,'String',size(h_T,1)*t_sample*1e3);
axes(handles.axes4);

% First shot 
set(handles.tval,'String',timei(ii)*1e3);
subplot(2,1,1), stem(tauv,abs(h_T(ii,:)))
xlabel('Delays [탎]')
ylabel('Amplitude')
grid on
subplot(2,1,2), stem(tauv,angle(h_T(ii,:)))
xlabel('Delays [탎]')
ylabel('Phase')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Outputs from this function are returned to the command line.
function varargout = Impulse_res_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global h_T;
global ii;
global timei;
global tauv;
set(handles.Start,'UserData',1);


while get(handles.Start,'Userdata')==1
    if ii<size(h_T,1)
        set(handles.tval,'String',timei(ii)*1e3);
        subplot(2,1,1), stem(tauv,abs(h_T(ii,:)))
        xlabel('Delays [탎]')
        ylabel('Amplitude')
        grid on
        subplot(2,1,2), stem(tauv,angle(h_T(ii,:)))
        xlabel('Delays [탎]')
        ylabel('Phase')
        grid on
        pause(0.1)
        ii=ii+1;
    else
        break
    end
end


% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
% hObject    handle to Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Start,'UserData',0);

% --- Executes on button press in Prev.
function Prev_Callback(hObject, eventdata, handles)
% hObject    handle to Prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global h_T;
global ii;
global timei;
global tauv;
set(handles.Start,'UserData',0);

if ii>1
    ii=ii-1;
    set(handles.tval,'String',timei(ii)*1e3);
    subplot(2,1,1), stem(tauv,abs(h_T(ii,:)))
    xlabel('Delays [탎]')
    ylabel('Amplitude')
    grid on
    subplot(2,1,2), stem(tauv,angle(h_T(ii,:)))
    xlabel('Delays [탎]')
    ylabel('Phase')
    grid on
else
    ii=1;
    set(handles.tval,'String',timei(ii)*1e3);
    subplot(2,1,1), stem(tauv,abs(h_T(ii,:)))
    xlabel('Delays [탎]')
    ylabel('Amplitude')
    grid on
    subplot(2,1,2), stem(tauv,angle(h_T(ii,:)))
    xlabel('Delays [탎]')
    ylabel('Phase')
    grid on
end

% --- Executes on button press in Next.
function Next_Callback(hObject, eventdata, handles)
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global h_T;
global ii;
global timei;
global tauv;
set(handles.Start,'UserData',0);
if ii<size(h_T,1)
    ii=ii+1;
    set(handles.tval,'String',timei(ii)*1e3);
    subplot(2,1,1), stem(tauv,abs(h_T(ii,:)))
    ylabel('Amplitude')
    xlabel('Delays [탎]')
    grid on
    subplot(2,1,2), stem(tauv,angle(h_T(ii,:)))
    ylabel('Phase')
    xlabel('Delays [탎]')
    grid on
end

% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global h_T;
global ii;
global timei;
global tauv;
set(handles.Start,'UserData',0);

ii=1;
set(handles.tval,'String',timei(ii)*1e3);
subplot(2,1,1), stem(tauv,abs(h_T(ii,:)))
xlabel('Delays [탎]')
ylabel('Amplitude')
grid on
subplot(2,1,2), stem(tauv,angle(h_T(ii,:)))
xlabel('Delays [탎]')
ylabel('Phase')
grid on


function tval_Callback(hObject, eventdata, handles)
% hObject    handle to tval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tval as text
%        str2double(get(hObject,'String')) returns contents of tval as a double


% --- Executes during object creation, after setting all properties.
function tval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fl_Callback(hObject, eventdata, handles)
% hObject    handle to fl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fl as text
%        str2double(get(hObject,'String')) returns contents of fl as a double


% --- Executes during object creation, after setting all properties.
function fl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_sample_Callback(hObject, eventdata, handles)
% hObject    handle to t_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_sample as text
%        str2double(get(hObject,'String')) returns contents of t_sample as a double


% --- Executes during object creation, after setting all properties.
function t_sample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
