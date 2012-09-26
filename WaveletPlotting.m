function varargout = WaveletPlotting(varargin)
% WAVELETPLOTTING MATLAB code for WaveletPlotting.fig
%      WAVELETPLOTTING, by itself, creates a new WAVELETPLOTTING or raises the existing
%      singleton*.
%
%      H = WAVELETPLOTTING returns the handle to a new WAVELETPLOTTING or the handle to
%      the existing singleton*.
%
%      WAVELETPLOTTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVELETPLOTTING.M with the given input arguments.
%
%      WAVELETPLOTTING('Property','Value',...) creates a new WAVELETPLOTTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WaveletPlotting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WaveletPlotting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WaveletPlotting

% Last Modified by GUIDE v2.5 25-Sep-2012 13:23:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WaveletPlotting_OpeningFcn, ...
    'gui_OutputFcn',  @WaveletPlotting_OutputFcn, ...
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


% --- Executes just before WaveletPlotting is made visible.
function WaveletPlotting_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WaveletPlotting (see VARARGIN)

% Choose default command line output for WaveletPlotting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Parse the varargin of the following format:
%value pairs (..., ..., 'Property', 'Value', ..., ...). There will be an
%even number.
for i=1:2:(length(varargin)-1)
    switch varargin{i}
        case 'cellTrackMatFilename'
            if exist(varargin{i+1},'file')
                handles.cellTrackMatFilename = varargin{i+1};
                load(handles.cellTrackMatFilename);
            end
        case 'signalMatFilename'
            if exist(varargin{i+1},'file')
                handles.signalMatFilename = varargin{i+1};
                load(handles.signalMatFilename);
                [pathstr, ~, ~] = fileparts(handles.signalMatFilename);
                cd(pathstr);
            end
        case 'mysignal'
            handles.signal = varargin{i+1};
        otherwise
    end
end

%%%%%% Assumptions
% 1. There is array with data from mysignal1
% 2. All cells have been sampled an equal number of times
% 3. There is no time vector with the time (days, hours, minutes) a sample was acquired.
% 4. The signal is positive, >= 0.

%%%%%% Variables
%handles.signal = mysignal1'; %a 2D array. 'size(handles.signal,1)' = # of cells. 'size(handles.signal,2)' = # of time points or samples.
numTime = size(handles.signal,2); % The number of time points in a signal
numCells = size(handles.signal,1); % The number of cells being analyzed
handles.signalInterpolated = cell(numCells,1); %see 1a
handles.signalPow2 = cell(numCells,1); %see 1b
handles.signalWindowed = cell(numCells,1); %see 1c
%fftLen
%timeDiff
%orignalSignalIndLeft
%percentSignalExpansion
%tukeyWindow
%halfSupport
%pow102
%wavelet_scales
%cwtftOutputTemp
%cwtftOutputTemp2
%handles.waveletScales
handles.waveletXfrm.type = 'mexh';
handles.waveletXfrm.cfs = cell(numCells,1);
handles.waveletXfrm.ridgemap = cell(numCells,1);
handles.waveletXfrm.ridgepeaks = cell(numCells,1);
handles.waveletXfrm(2).type = 'dog1';
handles.waveletXfrm(2).cfs = cell(numCells,1);
handles.waveletXfrm(2).ridgemap = cell(numCells,1);
handles.waveletXfrm(2).ridgepeaks = cell(numCells,1);
handles.waveletXfrm(3).type = 'morlex';
handles.waveletXfrm(3).cfs = cell(numCells,1);
handles.waveletXfrm(3).ridgemap = cell(numCells,1);
handles.waveletXfrm(3).ridgepeaks = cell(numCells,1);
handles.ind = 1;
handles.numCells = numCells;
%smoothsignal
%wavelet_peaks: the peaks found at each scale of a wavelet transform
%wavelet_peaks_alias: the labels assigned to each peak found in wavelet_peaks
%gap_limit: the number of scales that a peak in a ridge can skip over before breaking into a new ridge.
%ridge_map
penguinjet;
PCAin(numCells).method(3).low.num = [];
PCAin(numCells).method(3).low.meancfs = [];
PCAin(numCells).method(3).med.num = [];
PCAin(numCells).method(3).med.meancfs = [];
PCAin(numCells).method(3).high.num = [];
PCAin(numCells).method(3).high.meancfs = [];
PCAin(numCells).method(3).vector = [];
%meancfsall
%low
%med
%high
%ridgemapLogical
%ridgecfs
%%%%%% 1. pre-process data before performing the wavelet transform
% 1a. Interpolate between time points to fill in gaps; i.e. the sample must be uniformly sampled and if not must be made to appear so.
for i=1:numCells
    handles.signalInterpolated{i} = handles.signal(i,:);
end
% 1b. The wavelet transform can be computed using a FFT algorithm. In order to leverage the advantages of the FFT (speed), expand the number of data points to the nearest power of two of 150% of the original signal length. The 150% is needed for effective windowing.
fftLen= 2^(ceil(log(numTime*1.5)/log(2)));
timeDiff = fftLen - numTime;
for i=1:numCells
    if mod(timeDiff,2)
        %is odd
        handles.signalPow2{i} = padarray(handles.signalInterpolated{i},[0 (timeDiff-1)/2],'replicate');
        handles.signalPow2{i}(end+1) = handles.signalPow2{i}(end);
        originalSignalIndLeft = 1+(timeDiff-1)/2;
    else
        %is even
        handles.signalPow2{i} = padarray(handles.signalInterpolated{i},[0 timeDiff/2],'replicate');
        originalSignalIndLeft = 1+(timeDiff)/2;
    end
end
% 1c. Apply a Tukey window, a tapered cosine, to suppress edge effects. The expanded signal will be lowered to zero while the original signal remains the same.
percentSignalExpansion = timeDiff/fftLen;
tukeyWindow = tukeywin(fftLen, percentSignalExpansion);
for i=1:numCells
    handles.signalWindowed{i} = handles.signalPow2{i}.*tukeyWindow';
end

% 2. Calculate the wavelet transforms
%The cwtft function will only calculate coefficients for evenly spaced scales. However, this is an inconvenience when trying to look at trends that occur at different orders of magnitude. Therefore the cwtft function is called iteratively and then assembles the data into a struct that contains the coefficients and scales.

%Choosing the right scales to investigate can be a challenge, because it can feel subjective. One way is to include every integer scale up to the length of the signal, but this is probably too much information. Another is to choose scales on an exponential/log scale, but this might gloss over some important details. We'll use a hybrid between the two: filling in an exponential scale with uniform spacing. Hopefully this comprimise will deliver detail across several orders of magnitude.

%It was found to be that wavelet coefficients are no longer useful once the wavelet support is approx. half the length of the signal. The Mexican Hat wavelet has a support of 8 (or is it 11? check waveinfo('mexh')) at scale 1. The half scale is estimated to be the length of the signal/16 (assuming support of 8 time units).
% 2a. Determine the scales that will have wavelet coefficients.
halfSupport = numTime/16;
pow102 = ceil(log(halfSupport/10)/log(2));
if pow102<2
    pow102=2;
end
wavelet_scales = cell(1,pow102);
for i=1:pow102
    wavelet_scales{i} = (1:10)*2^(i-1)+10*(2^(i-1)-1);
end
if(wavelet_scales{pow102}(1)>halfSupport)
    wavelet_scales{pow102}(2:end) = [];
else
    wavelet_scales{pow102}(wavelet_scales{pow102}>halfSupport) = [];
end
handles.waveletScales = cell2mat(wavelet_scales);

% 2b. Calculate the mexican hat wavelet transform
for j=1:numCells
    cwtftOutputTemp = cwtft(handles.signalWindowed{j},'scales',wavelet_scales{1},'wavelet','mexh');
    cwtftOutputTemp2 = zeros(length(handles.waveletScales),fftLen);
    cwtftOutputTemp2(1:10,:) = cwtftOutputTemp.cfs;
    if (pow102>2)
        for i=2:(pow102-1)
            cwtftOutputTemp = cwtft(in,'scales',wavelet_scales{i},'wavelet','mexh');
            cwtftOutputTemp = cwtftOutputTemp.cfs;
            cwtftOutputTemp2((pow102-1)*10+1:(pow102-1)*10+10,:) = cwtftOutputTemp;
        end
    end
    cwtftOutputTemp = cwtft(handles.signalWindowed{j},'scales',wavelet_scales{pow102},'wavelet','mexh');
    cwtftOutputTemp = cwtftOutputTemp.cfs;
    cwtftOutputTemp2((end-length(wavelet_scales{pow102})+1):end,:) = cwtftOutputTemp;
    cwtftOutputTemp2 = real(cwtftOutputTemp2);
    for i = 1:length(handles.waveletScales)
        cwtftOutputTemp2(i,:) = cwtftOutputTemp2(i,:)/sqrt(handles.waveletScales(i));
    end
    handles.waveletXfrm(1).cfs{j} = cwtftOutputTemp2(:,originalSignalIndLeft:(originalSignalIndLeft+numTime-1));
end

% 2c. Calculate the derivative-of-gaussian wavelet transform
for j=1:numCells
    cwtftOutputTemp = cwtft(handles.signalWindowed{j},'scales',wavelet_scales{1},'wavelet',{'dog',1});
    cwtftOutputTemp2 = zeros(length(handles.waveletScales),fftLen);
    cwtftOutputTemp2(1:10,:) = cwtftOutputTemp.cfs;
    if (pow102>2)
        for i=2:(pow102-1)
            cwtftOutputTemp = cwtft(in,'scales',wavelet_scales{i},'wavelet',{'dog',1});
            cwtftOutputTemp = cwtftOutputTemp.cfs;
            cwtftOutputTemp2((pow102-1)*10+1:(pow102-1)*10+10,:) = cwtftOutputTemp;
        end
    end
    cwtftOutputTemp = cwtft(handles.signalWindowed{j},'scales',wavelet_scales{pow102},'wavelet',{'dog',1});
    cwtftOutputTemp = cwtftOutputTemp.cfs;
    cwtftOutputTemp2((end-length(wavelet_scales{pow102})+1):end,:) = cwtftOutputTemp;
    cwtftOutputTemp2 = real(cwtftOutputTemp2);
    for i = 1:length(handles.waveletScales)
        cwtftOutputTemp2(i,:) = cwtftOutputTemp2(i,:)/sqrt(handles.waveletScales(i));
    end
    handles.waveletXfrm(2).cfs{j} = cwtftOutputTemp2(:,originalSignalIndLeft:(originalSignalIndLeft+numTime-1));
end

% 2d. Calculate the non-analytic Morlet wavelet transform
for j=1:numCells
    cwtftOutputTemp = cwtft(handles.signalWindowed{j},'scales',wavelet_scales{1},'wavelet','morlex');
    cwtftOutputTemp2 = zeros(length(handles.waveletScales),fftLen);
    cwtftOutputTemp2(1:10,:) = cwtftOutputTemp.cfs;
    if (pow102>2)
        for i=2:(pow102-1)
            cwtftOutputTemp = cwtft(in,'scales',wavelet_scales{i},'wavelet','morlex');
            cwtftOutputTemp = cwtftOutputTemp.cfs;
            cwtftOutputTemp2((pow102-1)*10+1:(pow102-1)*10+10,:) = cwtftOutputTemp;
        end
    end
    cwtftOutputTemp = cwtft(handles.signalWindowed{j},'scales',wavelet_scales{pow102},'wavelet','morlex');
    cwtftOutputTemp = cwtftOutputTemp.cfs;
    cwtftOutputTemp2((end-length(wavelet_scales{pow102})+1):end,:) = cwtftOutputTemp;
    cwtftOutputTemp2 = real(cwtftOutputTemp2);
    for i = 1:length(handles.waveletScales)
        cwtftOutputTemp2(i,:) = cwtftOutputTemp2(i,:)/sqrt(handles.waveletScales(i));
    end
    handles.waveletXfrm(3).cfs{j} = cwtftOutputTemp2(:,originalSignalIndLeft:(originalSignalIndLeft+numTime-1));
end

% 3. Use ridge analysis to find peaks in the scalogram
% 3a. Find the 1% smoothend version of each signal.
smoothsignal = zeros(size(handles.signal));
for i = 1:numCells
    smoothsignal(i,:) = 0.01*smooth(handles.signal(i,:),10);
end
% 3b. Find the ridge map for each signal and wavelet transform
% 3b.i. Find the peaks within each scale
wavelet_peaks = cell(numCells,3);
for i = 1:numCells
    for j = 1:3
        wavelet_peaks{i,j} = cell(size(handles.waveletScales));
        for k = 1:length(handles.waveletScales)
            wavelet_peaks{i,j}{k} = first_pass_peak_detection(handles.waveletXfrm(j).cfs{i}(k,:), handles.waveletScales(k)*2+1, smoothsignal(i,:));
        end
    end
end
%3b.ii. Find the ridges
wavelet_peaks_alias = cell(numCells,3);
gap_limit = 2;
%%%%% ridge finding code
for y = 1:numCells
    for z = 1:3 %3 for the 3 types of wavelet transforms
        ridge_map = zeros(length(handles.waveletScales),numTime);
        wavelet_peaks_alias{y,z} = cell(size(handles.waveletScales));
        for i=1:length(wavelet_peaks_alias{y,z})
            wavelet_peaks_alias{y,z}{i} = zeros(size(wavelet_peaks{y,z}{i}));
        end
        for i=1:length(wavelet_peaks{y,z}{end})
            ridge_map(end,wavelet_peaks{y,z}{end}(i)) = i;
            wavelet_peaks_alias{y,z}{end}(i) = i;
        end
        ridge_counter = length(wavelet_peaks{y,z}{end}); %keeps track of the total number of ridges
        for i=size(handles.waveletXfrm(z).cfs{y},1):-1:(gap_limit+1)
            for j=1:length(wavelet_peaks{y,z}{i})
                for h=1:gap_limit
                    %Search for peaks within the window size for scale i.
                    low_bnd = wavelet_peaks{y,z}{i}(j) - 2*handles.waveletScales(i-h);
                    up_bnd = wavelet_peaks{y,z}{i}(j) + 2*handles.waveletScales(i-h);
                    low_set = wavelet_peaks{y,z}{i-h}>low_bnd;
                    up_set = wavelet_peaks{y,z}{i-h}<up_bnd;
                    %If a peak is found add it to the growing ridge
                    if any(low_set.*up_set)
                        ridge_set = low_set.*up_set;
                        if sum(ridge_set)==1
                            for k=1:length(ridge_set)
                                if ridge_set(k) && (wavelet_peaks_alias{y,z}{i-h}(k)==0 && ...
                                        sum(wavelet_peaks_alias{y,z}{i-h}==wavelet_peaks_alias{y,z}{i}(j))==0)
                                    ridge_map(i-h,wavelet_peaks{y,z}{i-h}(k)) = wavelet_peaks_alias{y,z}{i}(j);
                                    wavelet_peaks_alias{y,z}{i-h}(k) = wavelet_peaks_alias{y,z}{i}(j);
                                end
                            end
                        else
                            my_min = inf;
                            for k=1:length(ridge_set)
                                if ridge_set(k) && wavelet_peaks_alias{y,z}{i-h}(k)==0
                                    temp_ind=wavelet_peaks_alias{y,z}{i-(h-1)}==wavelet_peaks_alias{y,z}{i}(j);
                                    penultimate_ridge_position=wavelet_peaks{y,z}{i-(h-1)}(temp_ind);
                                    temp_min = abs(wavelet_peaks{y,z}{i-h}(k)-penultimate_ridge_position);
                                    if temp_min<my_min
                                        k_min = k;
                                        my_min = temp_min;
                                    end
                                end
                            end
                            if (sum(wavelet_peaks_alias{y,z}{i-h}==wavelet_peaks_alias{y,z}{i}(j))==0) && (k_min<=k)
                                ridge_map(i-h,wavelet_peaks{y,z}{i-h}(k_min)) = wavelet_peaks_alias{y,z}{i}(j);
                                wavelet_peaks_alias{y,z}{i-h}(k_min) = wavelet_peaks_alias{y,z}{i}(j);
                            end
                        end
                    end
                end
            end
            %Start new ridges for the peaks that were not assigned to ridges in the
            %scale below.
            new_ridge_set = wavelet_peaks_alias{y,z}{i-1}==0;
            for j=1:length(wavelet_peaks_alias{y,z}{i-1})
                if new_ridge_set(j)
                    ridge_counter = ridge_counter + 1;
                    wavelet_peaks_alias{y,z}{i-1}(j) = ridge_counter;
                    ridge_map(i-1,wavelet_peaks{y,z}{i-1}(j)) = ridge_counter;
                end
            end
        end
        handles.waveletXfrm(z).ridgemap{y} = ridge_map;
        low = 1:round(length(handles.waveletScales)/3);
        med = (low(end)+1):round(length(handles.waveletScales)*2/3);
        high = (med(end)+1):length(handles.waveletScales);
        ridgemapLogical = ridge_map;
        ridgemapLogical(ridgemapLogical>0)=1;
        PCAin(y).method(z).low.num = sum(sum(ridgemapLogical(low,:),2));
        PCAin(y).method(z).med.num = sum(sum(ridgemapLogical(med,:),2));
        PCAin(y).method(z).high.num = sum(sum(ridgemapLogical(high,:),2));
        ridgecfs = handles.waveletXfrm(z).cfs{y};
        ridgecfs(~ridgemapLogical) = 0;
        meancfsall = mean(ridgecfs,2);
        PCAin(y).method(z).low.meancfs = mean(meancfsall(low));
        PCAin(y).method(z).med.meancfs = mean(meancfsall(med));
        PCAin(y).method(z).high.meancfs = mean(meancfsall(high));
    end
    PCAin(y).vector = [PCAin(y).method(1).low.num, PCAin(y).method(1).low.meancfs, PCAin(y).method(1).med.num, PCAin(y).method(1).med.meancfs, PCAin(y).method(1).high.num, PCAin(y).method(1).high.meancfs, ...
        PCAin(y).method(2).low.num, PCAin(y).method(2).low.meancfs, PCAin(y).method(2).med.num, PCAin(y).method(2).med.meancfs, PCAin(y).method(2).high.num, PCAin(y).method(2).high.meancfs, ...
        PCAin(y).method(3).low.num, PCAin(y).method(3).low.meancfs, PCAin(y).method(3).med.num, PCAin(y).method(3).med.meancfs, PCAin(y).method(3).high.num, PCAin(y).method(3).high.meancfs];
end
save('PCAin.mat', 'PCAin');
PCAin2D = {PCAin(:).vector};
PCAin2D = cell2mat(PCAin2D');
save('PCAin2D.mat', 'PCAin2D');
%%%%% end of ridge finding code

% 4. Create features
% 4a. create features derived from the peaks identified on the ridgemap


guidata(hObject, handles);


% UIWAIT makes WaveletPlotting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WaveletPlotting_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in prev.
function prev_Callback(hObject, eventdata, handles)
% hObject    handle to prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ind~=1
    handles.ind = handles.ind - 1;
end
set(handles.edit_currentcell,'String',num2str(handles.ind));
wavelettype = get(handles.popupmenu_waveletoption,'Value');
cfsmap = handles.waveletXfrm(wavelettype).cfs{handles.ind};
cfsmap_min = min(min(cfsmap));
cfsmap = ((cfsmap - cfsmap_min)*253/(max(max(cfsmap))-cfsmap_min))+1;
cfsmap(handles.waveletXfrm(wavelettype).ridgemap{handles.ind}>0) = 255;
imagesc(cfsmap, 'Parent', handles.axes1);
axes(handles.axes2);
mysignal = handles.signal';
plot(mysignal(:,handles.ind));
guidata(hObject, handles);

% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.ind~=handles.numCells
    handles.ind = handles.ind + 1;
end
set(handles.edit_currentcell,'String',num2str(handles.ind));
wavelettype = get(handles.popupmenu_waveletoption,'Value');
cfsmap = handles.waveletXfrm(wavelettype).cfs{handles.ind};
cfsmap_min = min(min(cfsmap));
cfsmap = ((cfsmap - cfsmap_min)*253/(max(max(cfsmap))-cfsmap_min))+1;
cfsmap(handles.waveletXfrm(wavelettype).ridgemap{handles.ind}>0) = 255;
imagesc(cfsmap, 'Parent', handles.axes1);
axes(handles.axes2);
mysignal = handles.signal';
plot(mysignal(:,handles.ind));
guidata(hObject, handles);



function edit_currentcell_Callback(hObject, eventdata, handles)
% hObject    handle to edit_currentcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_currentcell as text
%        str2double(get(hObject,'String')) returns contents of edit_currentcell as a double
currentInd = str2num(get(handles.edit_currentcell,'String'));
handles.ind = currentInd;
wavelettype = get(handles.popupmenu_waveletoption,'Value');
cfsmap = handles.waveletXfrm(wavelettype).cfs{handles.ind};
cfsmap_min = min(min(cfsmap));
cfsmap = ((cfsmap - cfsmap_min)*253/(max(max(cfsmap))-cfsmap_min))+1;
cfsmap(handles.waveletXfrm(wavelettype).ridgemap{handles.ind}>0) = 255;
imagesc(cfsmap, 'Parent', handles.axes1);
axes(handles.axes2);
mysignal = handles.signal';
plot(mysignal(:,currentInd));
guidata(hObject, handles);








% --- Executes during object creation, after setting all properties.
function edit_currentcell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_currentcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_waveletoption.
function popupmenu_waveletoption_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_waveletoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_waveletoption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_waveletoption
currentInd = str2num(get(handles.edit_currentcell,'String'));
handles.ind = currentInd;
wavelettype = get(hObject,'Value');
imagesc(handles.waveletXfrm(wavelettype).cfs{handles.ind}, 'Parent', handles.axes1);
axes(handles.axes2);
mysignal = handles.signal';
assignin('base','mysignal',mysignal);
assignin('base','mywavelet',handles.waveletXfrm(wavelettype).cfs{handles.ind});
plot(mysignal(:,currentInd));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_waveletoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_waveletoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [out]=first_pass_peak_detection(x,winw_size,s)
%There are many ways to find peaks. The wavelet method is powerful at
%detecting peaks but is actually dependent on simpler peak detection
%methods.

%Simple Peak Finding Method: Scan a signal with a window. Find the max. If the max is greater
%than the left and right endpoints of the window it is a peak candidate.
%The window is centered at each point of a waveform, so a wider peak
%candidate will recieve more votes in a sense. If a window has more than
%one point with the max value then the left most index is used. A downside
%to this method is that it is sensitive to the size of the window. However,
%this seeming limitation is actually put to use in the wavelet method by
%scaling the window along with the wavelet.
length_x = length(x);
if size(x,2) == 1
    x=x';
end
if winw_size<3
    winw_size = 3;
    warning('win_size:small','The window size for first pass peak detection may be too small and therefore yield illogical results.');
end
pad=(winw_size-1)/2;
[x_padded,window_padded] = deal(zeros(1,length_x+winw_size-1));
window_padded(1:winw_size) = ones(1,winw_size);
x_padded(pad+1:end-pad) = x;
out=zeros(size(x));
for i=1:length_x
    seg = x_padded(logical(window_padded));
    [~,ind] = max(seg);
    if seg(ind) > seg(1) && seg(ind) > seg(end)
        out(i) = i + ind - pad - 1;
    end
    window_padded = circshift(window_padded,[0, 1]);
end
out(out==0) = [];
peak_candidates = unique(out);
peak_elected = peak_candidates; %This vector will be trimmed below
%Tally votes; if a peak candidate has less than or equal to pad votes they
%are disqualified
for i=length(peak_candidates):-1:1
    temp = out==peak_candidates(i);
    if sum(temp)<=pad/2
        peak_elected(i) = [];
    end
end
%Weed out unqualified peaks and peaks of questionable nature.
%If a peak has a negative wavelet coefficient value it is disqualified and
%if a peak is the first or last point of data it is disqualified.
for i=length(peak_elected):-1:1
    if x(peak_elected(i))<0
        peak_elected(i) = [];
    elseif peak_elected(i) == 1
        peak_elected(i) = [];
    elseif peak_elected(i) == length_x
        peak_elected(i) = [];
    end
end
%If a peak is less than the value of the 1% smoothened signal destroy its
%very existence. This is mainly to eliminate peaks that are very near to
%zero in the first few scales.
for i=length(peak_elected):-1:1
    if x(peak_elected(i))<s(peak_elected(i))
        peak_elected(i) = [];
    end
end
out = peak_elected;

function [] = penguinjet()
tuxedojet = [0,0,0;0,0,0.53125;0,0,0.546875;0,0,0.5625;0,0,0.578125;0,0,0.59375;0,0,0.609375;0,0,0.625;0,0,0.640625;0,0,0.65625;0,0,0.671875;0,0,0.6875;0,0,0.703125;0,0,0.71875;0,0,0.734375;0,0,0.75;0,0,0.765625;0,0,0.78125;0,0,0.796875;0,0,0.8125;0,0,0.828125;0,0,0.84375;0,0,0.859375;0,0,0.875;0,0,0.890625;0,0,0.90625;0,0,0.921875;0,0,0.9375;0,0,0.953125;0,0,0.96875;0,0,0.984375;0,0,1;0,0.015625,1;0,0.03125,1;0,0.046875,1;0,0.0625,1;0,0.078125,1;0,0.09375,1;0,0.109375,1;0,0.125,1;0,0.140625,1;0,0.15625,1;0,0.171875,1;0,0.1875,1;0,0.203125,1;0,0.21875,1;0,0.234375,1;0,0.25,1;0,0.265625,1;0,0.28125,1;0,0.296875,1;0,0.3125,1;0,0.328125,1;0,0.34375,1;0,0.359375,1;0,0.375,1;0,0.390625,1;0,0.40625,1;0,0.421875,1;0,0.4375,1;0,0.453125,1;0,0.46875,1;0,0.484375,1;0,0.5,1;0,0.515625,1;0,0.53125,1;0,0.546875,1;0,0.5625,1;0,0.578125,1;0,0.59375,1;0,0.609375,1;0,0.625,1;0,0.640625,1;0,0.65625,1;0,0.671875,1;0,0.6875,1;0,0.703125,1;0,0.71875,1;0,0.734375,1;0,0.75,1;0,0.765625,1;0,0.78125,1;0,0.796875,1;0,0.8125,1;0,0.828125,1;0,0.84375,1;0,0.859375,1;0,0.875,1;0,0.890625,1;0,0.90625,1;0,0.921875,1;0,0.9375,1;0,0.953125,1;0,0.96875,1;0,0.984375,1;0,1,1;0.015625,1,0.984375;0.03125,1,0.96875;0.046875,1,0.953125;0.0625,1,0.9375;0.078125,1,0.921875;0.09375,1,0.90625;0.109375,1,0.890625;0.125,1,0.875;0.140625,1,0.859375;0.15625,1,0.84375;0.171875,1,0.828125;0.1875,1,0.8125;0.203125,1,0.796875;0.21875,1,0.78125;0.234375,1,0.765625;0.25,1,0.75;0.265625,1,0.734375;0.28125,1,0.71875;0.296875,1,0.703125;0.3125,1,0.6875;0.328125,1,0.671875;0.34375,1,0.65625;0.359375,1,0.640625;0.375,1,0.625;0.390625,1,0.609375;0.40625,1,0.59375;0.421875,1,0.578125;0.4375,1,0.5625;0.453125,1,0.546875;0.46875,1,0.53125;0.484375,1,0.515625;0.5,1,0.5;0.515625,1,0.484375;0.53125,1,0.46875;0.546875,1,0.453125;0.5625,1,0.4375;0.578125,1,0.421875;0.59375,1,0.40625;0.609375,1,0.390625;0.625,1,0.375;0.640625,1,0.359375;0.65625,1,0.34375;0.671875,1,0.328125;0.6875,1,0.3125;0.703125,1,0.296875;0.71875,1,0.28125;0.734375,1,0.265625;0.75,1,0.25;0.765625,1,0.234375;0.78125,1,0.21875;0.796875,1,0.203125;0.8125,1,0.1875;0.828125,1,0.171875;0.84375,1,0.15625;0.859375,1,0.140625;0.875,1,0.125;0.890625,1,0.109375;0.90625,1,0.09375;0.921875,1,0.078125;0.9375,1,0.0625;0.953125,1,0.046875;0.96875,1,0.03125;0.984375,1,0.015625;1,1,0;1,0.984375,0;1,0.96875,0;1,0.953125,0;1,0.9375,0;1,0.921875,0;1,0.90625,0;1,0.890625,0;1,0.875,0;1,0.859375,0;1,0.84375,0;1,0.828125,0;1,0.8125,0;1,0.796875,0;1,0.78125,0;1,0.765625,0;1,0.75,0;1,0.734375,0;1,0.71875,0;1,0.703125,0;1,0.6875,0;1,0.671875,0;1,0.65625,0;1,0.640625,0;1,0.625,0;1,0.609375,0;1,0.59375,0;1,0.578125,0;1,0.5625,0;1,0.546875,0;1,0.53125,0;1,0.515625,0;1,0.5,0;1,0.484375,0;1,0.46875,0;1,0.453125,0;1,0.4375,0;1,0.421875,0;1,0.40625,0;1,0.390625,0;1,0.375,0;1,0.359375,0;1,0.34375,0;1,0.328125,0;1,0.3125,0;1,0.296875,0;1,0.28125,0;1,0.265625,0;1,0.25,0;1,0.234375,0;1,0.21875,0;1,0.203125,0;1,0.1875,0;1,0.171875,0;1,0.15625,0;1,0.140625,0;1,0.125,0;1,0.109375,0;1,0.09375,0;1,0.078125,0;1,0.0625,0;1,0.046875,0;1,0.03125,0;1,0.015625,0;1,0,0;0.984375,0,0;0.96875,0,0;0.953125,0,0;0.9375,0,0;0.921875,0,0;0.90625,0,0;0.890625,0,0;0.875,0,0;0.859375,0,0;0.84375,0,0;0.828125,0,0;0.8125,0,0;0.796875,0,0;0.78125,0,0;0.765625,0,0;0.75,0,0;0.734375,0,0;0.71875,0,0;0.703125,0,0;0.6875,0,0;0.671875,0,0;0.65625,0,0;0.640625,0,0;0.625,0,0;0.609375,0,0;0.59375,0,0;0.578125,0,0;0.5625,0,0;0.546875,0,0;0.53125,0,0;0.515625,0,0;1,1,1];
colormap(tuxedojet);
