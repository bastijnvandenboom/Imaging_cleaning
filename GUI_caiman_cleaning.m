% GUI to plot caiman cnmf output and manual curate ROIs

% GUI has 2 main functions, delete ROIs and merge ROIs

% caiman/cnmf rejects many ROI (rois_bad), you can decide to plot those
% (check plot bad rois). if you do, those ROIs will go to the delete list

% due to it's underlying structure, you need to delete or merge (not both
% at the same time)

% you can perform a delete after a merge and the otherway around, except
% that merge pairs will not be redetected -> reload data to do so

% in practise, you might have multiple rounds of deleting and merging
% (run GUI, save data, restart GUI, reload those data, and continue)

% based on miniscope GUI (Aishwarya Parthasarathy and Bastijn van den Boom
% at Willuhn lab, NIN)
% 2025 Bastijn van den Boom (Sabatini lab, HMS)


function varargout = GUI_caiman_cleaning(varargin)
% GUI_CAIMAN_CLEANING MATLAB code for GUI_caiman_cleaning.fig
%      GUI_CAIMAN_CLEANING, by itself, creates a new GUI_CAIMAN_CLEANING or raises the existing
%      singleton*.
%
%      H = GUI_CAIMAN_CLEANING returns the handle to a new GUI_CAIMAN_CLEANING or the handle to
%      the existing singleton*.
%
%      GUI_CAIMAN_CLEANING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CAIMAN_CLEANING.M with the give0n input arguments.
%
%      GUI_CAIMAN_CLEANING('Property','Value',...) creates a new GUI_CAIMAN_CLEANING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_caiman_cleaning_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_caiman_cleaning_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_caiman_cleaning

% Last Modified by GUIDE v2.5 07-Feb-2025 15:40:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_caiman_cleaning_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_caiman_cleaning_OutputFcn, ...
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


% --- Executes just before GUI_caiman_cleaning is made visible.
function GUI_caiman_cleaning_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_caiman_cleaning (see VARARGIN)

% Choose default command line output for GUI_caiman_cleaning
handles.output = hObject;

% allow updating of plotting bad ROIs
handles.plot_bad_update = 1;

% allow zooming figures
zoom on

% start_gui handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_caiman_cleaning wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_caiman_cleaning_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BASIC GUI


% --- Executes on button press in LoadFile.
function LoadFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clear command window
clc

% clear handles fields of caiman
handles.rawdata1 = [];
handles.updatedCraw = [];
handles.updatedC = [];
handles.updatedA = [];
handles.updatedS = [];
handles.updatedresiduals = [];
handles.updatedSNRcnmf = [];
handles.SNR = [];
handles.mergecells = [];
handles.delcells = [];
handles.pixels_A = [];
handles.dist_mat = [];
handles.corr_mat = [];
handles.dist = [];
handles.corr = [];
handles.pairs_currentcompare = [];
handles.multimergeidx = [];
handles.plot_roi = [];
handles.manualpair = [];
set(handles.cell_pair1, 'Value', 0);
set(handles.cell_pair2, 'Value', 0);
set(handles.cell_pair_dist, 'Value', 0);
set(handles.cell_pair_corr, 'Value', 0);
set(handles.n_pair, 'String', []);
set(handles.plot_rois, 'Value', 0);
set(handles.roi_idx, 'Value', 0);
set(handles.n_rois, 'Value', 0);
set(handles.roi_snr, 'Value', 0);
set(handles.roi_size, 'Value', 0);

% clear panels
set(handles.deletelist, 'String', sprintf('ROIs to delete'), 'Value', 1); % update values in deletelist
set(handles.mergelist, 'String', sprintf('ROIs to merge'), 'Value', 1); % update values in mergelist
set(handles.checkmergetxt, 'String', 'Click "Find repeats"'); % update values in mergelist

% user select file
[filename1, filepath1] = uigetfile('*.*', 'Select raw .hdf5 caiman (cnmf) output file or preprocessed .mat file');
cd(filepath1)

% check if we want to plot bad ROIs
handles.plot_bad = get(handles.check_plot_bad, 'Value'); % 0=good ROIs, 1=bad ROIs
handles.plot_bad_update = 0;    % prevent updating bad ROI plotting

% check if file is a (preprocessed) .mat file or directly .hdf5 output of caiman (cnmf)
tmp = strsplit(filename1, '.');
tmp = tmp{end};

% update user and load data
if strcmp(tmp, 'mat') % mat file
    % update user
    set(handles.user_alert, 'String', sprintf('Loading .mat file, please wait'));
    % load file
    handles.rawdata1 = load([filepath1 filename1]);   % all raw data
    % update bad ROIs if needed
    if handles.plot_bad == 1 % update results
        handles.rawdata1.results.C_raw = handles.rawdata1.results.raw.C_raw;
        handles.rawdata1.results.C = handles.rawdata1.results.raw.C;
        handles.rawdata1.results.S = handles.rawdata1.results.raw.S;
        handles.rawdata1.results.A = handles.rawdata1.results.raw.A;
        handles.rawdata1.results.rois_bad = handles.rawdata1.results.raw.rois_bad;
        try % not everyone has these vars
            handles.rawdata1.results.residuals = handles.rawdata1.results.raw.residuals;
            handles.rawdata1.results.SNR_cnmf = handles.rawdata1.results.raw.SNR_cnmf;
        catch end
    end

elseif strcmp(tmp, 'hdf5') % raw output
    % update user
    set(handles.user_alert, 'String', sprintf('Loading .hdf5 file, please wait'));
    % load file
    handles.rawdata1.results = load_cnmf_hdf5(filename1, filepath1, handles);

else
    % update user
    set(handles.user_alert, 'String', sprintf('Unknown file format, please select .mat or .hdf5'));
    error('Unknown file format')
end

% extract variables
handles.updatedCraw = handles.rawdata1.results.C_raw;
handles.updatedC = handles.rawdata1.results.C;
handles.updatedA = handles.rawdata1.results.A;
handles.updatedS = handles.rawdata1.results.S;
try % not everyone has these vars 
    handles.updatedresiduals = handles.rawdata1.results.residuals;
    handles.updatedSNRcnmf = handles.rawdata1.results.SNR_cnmf;
catch end
try % check if we've ran the GUI before, important for rois_bad!
    handles.GUI = handles.rawdata1.results.GUI;
catch
    handles.GUI = 'false';
end

% check if we need to update bad ROIs
if handles.plot_bad == 1 % update results
    if strcmp(handles.GUI, 'true') % have we ran the GUI before?
        % update user
        set(handles.user_alert, 'String', sprintf('GUI has been used before on bad rois!'));
        error('GUI has been used before on these data, you cannot use bad ROIs because the ROI numbers probably will not match')
    end
    handles.rois_bad = handles.rawdata1.results.raw.rois_bad;
end

% update total n ROIs
set(handles.n_rois, 'String', sprintf('%d', size(handles.updatedCraw,1))); % update total n ROIs

% update user
set(handles.user_alert, 'String', sprintf('Data loaded, start GUI'));

% store data
guidata(hObject,handles);


% --- Executes on button press in start_gui.
function start_gui_Callback(hObject, eventdata, handles)
% hObject    handle to start_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update user
set(handles.user_alert, 'String', sprintf('Starting GUI, please wait'));
pause(0.1); % make sure to update user_alert

% clear panels
set(handles.deletelist, 'String', sprintf('ROIs to delete'), 'Value', 1); % update values in deletelist
set(handles.mergelist, 'String', sprintf('ROIs to merge'), 'Value', 1); % update values in mergelist
set(handles.checkmergetxt, 'String', 'Click "Find repeats"'); % update values in mergelist

% move all the bad ROIs to the delete list if we're plotting bad ROIs
if handles.plot_bad == 1
    handles.rois_bad = double(handles.rois_bad); % make it double
    handles.delcells = handles.rois_bad; % update delcells list
    set(handles.deletelist, 'String', handles.rois_bad, 'Value', length(handles.rois_bad)); % add values to delete list
end

% define colors
handles.colors = [56/255 176/255 0/255; 0/255 0/255 255/255; 255/255 0/255 0/255]; % green, blue, red [single ROI, 2 ROIs]
handles.color_plots = [.5 0 .5; 251/255 111/255 146/255]; % purple, pink [distance x corr and histo plots, selected ROIs]

% calculate SNR
handles.SNR = get_SNR(handles.updatedC, handles.updatedCraw);

% % find peaks and baseline
% for i = 1:size(handles.rawdata1.results.C_raw,1)    % per roi
%     [pk{i}, loc{i}] = findpeaks(handles.updatedC(i,:)); % , 3
%
% end
% % find baseline
%
% % calculate SNR - divide by baseline
% handles.SNR=[];deletelist
% for i=1:size(handles.rawdata1.results.C_raw,1)    % per roi [size(handles.updatedCraw,1)]
%     [b, sn] = estimate_baseline_noise(handles.updatedCraw(i,:))
%     handles.SNR(i) = median(pk{i})/sn;
%     clear sn;
% end

% % find ROIs with low SNR and remove those ROIshandles.SNR
% idx_snr = find(handles.SNR < str2num(get(handles.snr_thres,'String')));
% handles.updatedCraw(idx_snr,:) = [];
% handles.updatedC(idx_snr,:) = [];
% handles.updatedS(idx_snr,:) = [];
% handles.updatedA(:,idx_snr) = [];
% handles.SNR(idx_snr) = [];

% find ROIs with low SNR and add to deletelist
idx_snr = find(handles.SNR < str2num(get(handles.snr_thres,'String')));
if ~isempty(idx_snr)
    tmp_list = handles.delcells; % get list
    tmp_list = [tmp_list idx_snr]; % add roi
    handles.delcells = tmp_list; % update del list
    set(handles.deletelist, 'String', handles.delcells, 'Value', length(handles.delcells)); % update values in list
end

% get full ROI, center of mass, and size A
[Cent_N, handles.pixels_A, handles.Afull] = get_shape_size_A(handles);

% loop through roi pairs to get correlation and distance - SLOW step
for i = 1:size(handles.updatedCraw,1) % for each roi
    for j = i:size(handles.updatedCraw,1) % for each other roi
        handles.dist_mat(j,i) = sqrt((Cent_N(i,1) - Cent_N(j,1))^2 + (Cent_N(i,2) - Cent_N(j,2))^2); % get distance
        handles.corr_mat(j,i) = corr(handles.updatedCraw(i,:)',handles.updatedCraw(j,:)'); % get correlation
    end
end

% for dist and corr: get lower triangular part of matrix by replacing the top part with NaNs
tril_mask = tril(true(size(handles.dist_mat)), -1);
handles.dist_mat(~tril_mask) = NaN;
handles.corr_mat(~tril_mask) = NaN;

% find ROIs that are too small (a_min) or big (a_max)
idx_amin = [];
idx_amax = [];
if str2num(get(handles.a_min,'String')) ~= 0
    idx_amin = find(handles.pixels_A < str2num(get(handles.a_min,'String')));
end
if str2num(get(handles.a_max,'String')) ~= 0
    idx_amax = find(handles.pixels_A > str2num(get(handles.a_max,'String')));
end
idx_a = unique([idx_amin idx_amax]);

% % remove those ROIs
% handles.updatedCraw(idx_a,:) = [];
% handles.updatedC(idx_a,:) = [];
% handles.updatedS(idx_a,:) = [];
% handles.updatedA(:,idx_a) = [];
% handles.SNR(idx_a) = [];

% add those ROIs to the dellist
if ~isempty(idx_a)
    tmp_list = handles.delcells; % get list
    tmp_list = [tmp_list idx_a]; % add roi
    handles.delcells = tmp_list; % update del list
    set(handles.deletelist, 'String', handles.delcells, 'Value', length(handles.delcells)); % update values in list
end

% find pairs of ROIs that fullful dist_thres & corr_thres
[handles.pairs_cells1, handles.pairs_cells2] = find(handles.dist_mat < str2num(get(handles.dist_thres,'String')) & ...
    handles.corr_mat > str2num(get(handles.corr_thres,'String')));

% set dropdown to merge pairs
% set(handles.roi1, 'String', unique([handles.pairs_cells1 handles.pairs_cells2]));
% set(handles.roi2, 'String', unique([handles.pairs_cells1 handles.pairs_cells2]));

% set dropdown to all possible rois
set(handles.roi1, 'String', [1:size(handles.updatedCraw,1)]);
set(handles.roi2, 'String', [1:size(handles.updatedCraw,1)]);

% make matrix into vector
handles.dist = reshape(handles.dist_mat, 1, size(handles.dist_mat,1)*size(handles.dist_mat,2));
handles.corr = reshape(handles.corr_mat, 1, size(handles.corr_mat,1)*size(handles.corr_mat,2));
idx = isnan(handles.dist); % find NaN's (which is top part and the diagonal)
handles.dist(idx)=[];
handles.corr(idx)=[];

% if we include bad ROIs, we'll start by plotting single ROIs. otherwise ROI pairs
if handles.plot_bad == 1 % bad ROIs
    % store that we plot a single ROI
    handles.to_plot = 1;

    % ROI to plot
    tmp = str2num(get(handles.deletelist,'String'));
    idx = tmp(1);
    handles.plot_roi = idx;

    % select the first ROI in the dellist
    set(handles.deletelist, 'String', get(handles.deletelist,'String'), 'Value', 1);

    % start a compare pair
    handles.pairs_currentcompare = 1;

    % update user
    set(handles.user_alert, 'String', sprintf('Delete or merge, then save'));
else
    % store that we plot a mergepair
    handles.to_plot = 0;

    % ROI(s) to plot
    handles.pairs_currentcompare = 1;

    % based on the settings, we might find no pairs at all
    try
        idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2
    catch
        % store that we plot a single ROI since there are no merge pairs
        handles.to_plot = 1;
        % plot the first ROI
        idx = 1;
        % update user
        set(handles.user_alert, 'String', sprintf('No merge pairs found!\nDelete or merge, then save'));
    end
end

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% define slider location
set(handles.slidertoaddtodellist, 'Min', 1, 'Max', size(handles.updatedCraw,1), 'SliderStep', [1/(size(handles.updatedCraw,1)-1),1], 'Value', 1)

% plot number of roi pairs that fulfill criteria
set(handles.n_pair, 'String', sprintf('%d pairs of ROIs satisfy the thresholds', length(handles.pairs_cells1)));

% plot histogram of SNR
plot_histo_snr(handles)

% plot histogram of size
plot_histo_size(handles)

% store data
guidata(hObject,handles);


% --- Executes on button press in restart_gui.
function restart_gui_Callback(hObject, eventdata, handles)
% hObject    handle to restart_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current figure handle
fig = handles.figure1;  % Adjust if your figure has a different name

% Get the GUI's filename (without the .fig extension)
guiName = 'GUI_caiman_cleaning';  % Change this to your actual GUI function name

% Close the current GUI
close(fig);

% Restart the GUI
feval(guiName);


% --- Executes on button press in del.
function del_Callback(hObject, eventdata, handles)
% hObject    handle to del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get deletelist and delete values
del_rois = str2num(get(handles.deletelist,'String'));
handles.updatedCraw(del_rois,:) = [];
handles.updatedC(del_rois,:) = [];
handles.updatedS(del_rois,:) = [];
handles.updatedA(:,del_rois) = [];
try % not everyone has these vars
    handles.updatedresiduals(del_rois,:) = [];
    handles.updatedSNRcnmf(del_rois,:) = [];
catch end
handles.SNR(:,del_rois) = [];
handles.pixels_A(:,del_rois) = [];

% clear del_list
set(handles.deletelist,'String',[]);

% update number of ROIs
set(handles.user_alert, 'String', sprintf('Done deleting, save data!\n Number of ROIs: %d', size(handles.updatedCraw,1)))

% update total n ROIs
set(handles.n_rois, 'String', sprintf('%d', size(handles.updatedCraw,1))); % update total n ROIs

% delete rois_bad values if we plotted those
if handles.plot_bad == 1 % plot bad ROIs
    handles.rois_bad = [];
end

% get full ROI, center of mass, and size A
[~, handles.pixels_A, handles.Afull] = get_shape_size_A(handles);

% plot histogram of SNR
plot_histo_snr(handles)

% plot histogram of size
plot_histo_size(handles)

% store data
guidata(hObject,handles);


% --- Executes on button press in merge.
function merge_Callback(hObject, eventdata, handles)
% hObject    handle to merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get ROIs to merge
merge_rois = get(handles.mergelist,'String');

% average activity (Craw, C, S) and combine contours (A)
rois_keep = []; % first ROIs of the pairs, to  keep
rois_delete = []; % all the other ROIs, to delete
for i = 1:length(merge_rois) % per mergepair

    % temporal activity
    rois_pair = sort(str2num(cell2mat(merge_rois(i))),'ascend');
    rois_craw = [];
    rois_c = [];
    rois_s = [];
    rois_res = [];
    rois_SNRcnmf = [];
    for u = 1:length(rois_pair) % per roi in mergepair
        rois_craw = [rois_craw; handles.updatedCraw(rois_pair(u),:)];
        rois_c = [rois_c; handles.updatedC(rois_pair(u),:)];
        rois_s = [rois_s; handles.updatedS(rois_pair(u),:)];
        try % not everyone has these vars
            rois_res = [rois_res; handles.updatedresiduals(rois_pair(u),:)];
            rois_SNRcnmf = [rois_SNRcnmf; handles.updatedSNRcnmf(rois_pair(u),:)];
        catch end
    end
    merged_craw(i,:) = mean(rois_craw);
    merged_c(i,:) = mean(rois_c);
    merged_s(i,:) = mean(rois_s);
    merged_res(i,:) = mean(rois_res);
    merged_SNRcnmf(i,:) = mean(rois_SNRcnmf);

    % spatial contour
    rois_a = [];
    % merged_contour(:,i) = handles.updatedA(:,rois_pair(1));
    for u = 1:length(rois_pair) % per roi in mergepair
        rois_a = [rois_a handles.updatedA(:,rois_pair(u))];
    end
    merged_a(:,i) = sum(rois_a, 2);

    % keep track of ROIs
    rois_keep = [rois_keep rois_pair(1)]; % list of all first ROIs
    rois_delete = [rois_delete rois_pair(2:end)]; % list of all second and more ROIs
end
handles.updatedCraw(rois_keep,:) = merged_craw; % store data on position first ROIs
handles.updatedCraw(rois_delete,:) = []; % remove other ROIs
handles.updatedC(rois_keep,:) = merged_c; % store data on position first ROIs
handles.updatedC(rois_delete,:) = []; % remove other ROIs
handles.updatedS(rois_keep,:) = merged_s; % store data on position first ROIs
handles.updatedS(rois_delete,:) = []; % remove other ROIs
handles.updatedA(:,rois_keep) = merged_a; % same for A
handles.updatedA(:,rois_delete) = [];
try % not everyone has these vars
    handles.updatedresiduals(rois_keep,:) = merged_res; % same for residuals
    handles.updatedresiduals(rois_delete,:) = [];
    handles.updatedSNRcnmf(rois_keep,:) = merged_SNRcnmf; % same for SNR_cnmf
    handles.updatedSNRcnmf(rois_delete,:) = [];
catch end

% clear merge data
set(handles.mergelist,'String',[]);
handles.pwdist = [];
handles.crosscoef = [];
handles.pairs_cells1 = [];
handles.pairs_cells2 = [];
handles.mergecells = [];

% update user
set(handles.user_alert, 'String', sprintf('Done merging, save data!\n Number of ROIs: %d', size(handles.updatedCraw,1)))

% update total n ROIs
set(handles.n_rois, 'String', sprintf('%d', size(handles.updatedCraw,1))); % update total n ROIs

% delete rois_bad values if we plotted those
if handles.plot_bad == 1 % plot bad ROIs
    handles.rois_bad = [];
end

% get full ROI, center of mass, and size A
[~, handles.pixels_A, handles.Afull] = get_shape_size_A(handles);

% plot histogram of SNR
plot_histo_snr(handles)

% plot histogram of size
plot_histo_size(handles)

% store data
guidata(hObject,handles);


function user_alert_Callback(hObject, eventdata, handles)
% hObject    handle to user_alert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of user_alert as text
%        str2double(get(hObject,'String')) returns contents of user_alert as a double


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% calculate SNR on the left over neurons
handles.SNR = get_SNR(handles.updatedC, handles.updatedCraw);

% find folder to save and make name
folder_name = uigetdir;
file_str = sprintf('%s\\updated_caiman.mat', folder_name);
  
% update user
set(handles.user_alert, 'String', sprintf('Storing data, please wait'));

% generate results struct
results = handles.rawdata1.results;
results.C_raw = handles.updatedCraw;
results.C = handles.updatedC;
results.S = handles.updatedS;
results.A = handles.updatedA;
try % not everyone has these vars
    results.residuals = handles.updatedresiduals;
    results.SNR_cnmf = handles.updatedSNRcnmf;
catch end
results.SNR_gui = handles.SNR';
results.GUI = 'true';

% allow loading bad ROIs
handles.plot_bad_update = 1;    % allow updating bad ROI plotting

% ROI(s) to plot
idx = 1:size(results.A,2); % all ROIs

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% save data
save(file_str,'results');

% save screenshot figure
fig1_str = sprintf('%s\\screenshot_GUI.png', folder_name);
frame = getframe(handles.figure1); % Capture the GUI figure
imwrite(frame.cdata, fig1_str); % Save as an image file

% feedback to user
set(handles.user_alert, 'String', sprintf('Done saving! Saved at: %s', file_str));

% go up one folder (store data in new folder)
% cd ..

% store data
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELETE ROIs


% --- Executes on slider movement.
function slidertoaddtodellist_Callback(hObject, eventdata, handles)
% hObject    handle to slidertoaddtodellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% store that we plot 1 roi
handles.to_plot = 1;

% define range of slider
set(handles.slidertoaddtodellist, 'Max', size(handles.updatedCraw,1), 'Value', round(get(handles.slidertoaddtodellist,'Value')), 'SliderStep', [1/(size(handles.updatedCraw,1)-1),1]);

% ROI(s) to plot
handles.plot_roi = get(handles.slidertoaddtodellist,'Value');
idx = handles.plot_roi;

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on button press in add2dellist.
function add2dellist_Callback(hObject, eventdata, handles)
% hObject    handle to add2dellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store that we plot 1 roi
handles.to_plot = 1;

% add to del list
tmp_list = handles.delcells; % get list
tmp_list = [tmp_list handles.plot_roi]; % add roi
handles.delcells = tmp_list; % update del list
set(handles.deletelist, 'String', tmp_list, 'Value', length(tmp_list)); % update values in list

% add to go to next ROI (unless it's the last ROI)
if handles.plot_roi ~= size(handles.updatedCraw,1)

    % ROI(s) to plot
    handles.plot_roi= handles.plot_roi + 1;
    idx = handles.plot_roi;
    set(handles.slidertoaddtodellist, 'Max', size(handles.updatedCraw,1), 'Value', handles.plot_roi, 'SliderStep', [1/(size(handles.updatedCraw,1)-1),1]);

    % plot background (Cn or Mn) and ROI contours (idx)
    plot_spatial_components(handles, idx)

    % plot temporal traces of the first pair
    plot_temporal_traces(handles, idx)

    % plot ROI pair info
    plot_roi_info(handles, idx)

end

% store data
guidata(hObject,handles);


% --- Executes on button press in rem4dellist.
function rem4dellist_Callback(hObject, eventdata, handles)
% hObject    handle to rem4dellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store that we plot 1 roi
handles.to_plot = 1;

% remove roi from dellist
k1 = get(handles.deletelist,'Value'); % get location
k2 = get(handles.deletelist,'String'); % get total list
k2(k1,:) = []; % remove from list
handles.delcells(k1) = []; % update delcells list
set(handles.deletelist, 'String', k2, 'Value', size(k2,1)); % update deletelist

% plot next deletelist ROI (unless it's the last ROI)
if k1 <= size(k2,1) % select the next ROI

    % ROI(s) to plot
    idx = str2num(k2(k1,:));
    handles.plot_roi = idx;

    % select the next ROI in the dellist
    set(handles.deletelist, 'String', k2, 'Value', k1);

elseif size(k2,1) == 0 % we're removed the last ROI in the deletelist

    % ROI(s) to plot -> make it ROI 1
    idx = 1;
    handles.plot_roi = idx;

elseif k1-1 == size(k2,1) % we removed the last ROI

    % ROI(s) to plot
    idx = str2num(k2(k1-1,:));
    handles.plot_roi = idx;

    % select the next ROI in the dellist
    set(handles.deletelist, 'String', k2, 'Value', k1-1);

end

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on selection change in deletelist.
function deletelist_Callback(hObject, eventdata, handles)
% hObject    handle to deletelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns deletelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from deletelist

% store that we plot 1 roi
handles.to_plot = 1;

% ROI(s) to plot
% find value of roi in deletelist
k1 = get(handles.deletelist,'Value');
k2 = str2num(get(handles.deletelist,'String'));
idx = k2(k1);
handles.plot_roi = idx;

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% store data
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MERGE ROI PAIRS

function n_pair_Callback(hObject, eventdata, handles)
% hObject    handle to n_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_pair as text
%        str2double(get(hObject,'String')) returns contents of n_pair as a double


% --- Executes on selection change in roi1.
function roi1_Callback(hObject, eventdata, handles)
% hObject    handle to roi1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roi1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roi1

% % assign all ROIs to roi1
% set(handles.roi1, 'String', unique([handles.pairs_cells1 handles.pairs_cells2]));

% get current rois
if handles.to_plot ~= 2 % we're not trying to manually plot, yet
    try % we might plot a single ROI
        idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2
    catch
        idx = [1 handles.plot_roi]; % just set ROI 1 on position 1
    end
else
    idx = [handles.manualpair]; % get manual ROI 1 and ROI 2 values
end

% store that we plot a manual pair
handles.to_plot = 2;

% replace one with new roi
roi_1_val = get(handles.roi1, 'Value'); % get roi1 value
idx = [roi_1_val idx(2)];
handles.manualpair = idx;

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on selection change in roi2.
function roi2_Callback(hObject, eventdata, handles)
% hObject    handle to roi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roi2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roi2

% % assign all ROIs to roi2
% set(handles.roi2,'String',unique([handles.pairs_cells1 handles.pairs_cells2]));

% get current rois
if handles.to_plot ~= 2 % we're not trying to manually plot, yet
    try % we might plot a single ROI
        idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2
    catch
        idx = [handles.plot_roi 1]; % just set ROI 1 on position 2
    end
else
    idx = [handles.manualpair]; % get manual ROI 1 and ROI 2 values
end

% store that we plot a manual pair
handles.to_plot = 2;

% replace one with new roi
roi_2_val = get(handles.roi2, 'Value'); % get roi2 value
idx = [idx(1) roi_2_val];
handles.manualpair = idx;

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on button press in add2merge.
function add2merge_Callback(hObject, eventdata, handles)
% hObject    handle to add2merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store that we plot a mergepair
handles.to_plot = 0;

% ROIs to store and plot - only pairs
% un_cells = unique([handles.pairs_cells1 handles.pairs_cells2]); % get unique pairs (cells1 and cells2 should be the same though)
% roi_1_val = un_cells(get(handles.roi1, 'Value')); % get roi1 value
% roi_2_val = un_cells(get(handles.roi2, 'Value')); % get roi2 value

% ROIs to store and plot
roi_vals = get(handles.roi1, 'Value'); % get roi1 value
roi_vals = sort([roi_vals get(handles.roi2, 'Value')], 'descend'); % get roi2 value

% update mergelist
mlist = handles.mergecells; % get mergelist
mlist = [mlist; {mat2str(roi_vals)} ]; % make updated list
set(handles.mergelist, 'String', mlist, 'Value', length(mlist)) % update list GUI
handles.mergecells = mlist; % update actual list

% ROI(s) to plot
idx = roi_vals; % ROI 1 and ROI 2

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on button press in pre.
function pre_Callback(hObject, eventdata, handles)
% hObject    handle to pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store that we plot a mergepair
handles.to_plot = 0;

% go to the previous pair (if this was not the first pair)
if handles.pairs_currentcompare~=1

    % roi pair to plot
    handles.pairs_currentcompare = handles.pairs_currentcompare - 1;

    % ROI(s) to plot
    idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2

    % plot background (Cn or Mn) and ROI contours (idx)
    plot_spatial_components(handles, idx)

    % plot temporal traces of the first pair
    plot_temporal_traces(handles, idx)

    % plot ROI pair info
    plot_roi_info(handles, idx)
end

% store data
guidata(hObject,handles);


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store that we plot a mergepair
handles.to_plot = 0;

% go to the next pair (if this was not the last pair)
if handles.pairs_currentcompare~=length(handles.pairs_cells1)

    % roi pair to plot
    handles.pairs_currentcompare = handles.pairs_currentcompare + 1;

    % ROI(s) to plot
    idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2

    % plot background (Cn or Mn) and ROI contours (idx)
    plot_spatial_components(handles, idx)

    % plot temporal traces of the first pair
    plot_temporal_traces(handles, idx)

    % plot ROI pair info
    plot_roi_info(handles, idx)
end

% store data
guidata(hObject,handles);


% --- Executes on button press in addlazy.
function addlazy_Callback(hObject, eventdata, handles)
% hObject    handle to addlazy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store that we plot a mergepair
handles.to_plot = 0;

% add current pair to mergelist
tt = handles.pairs_currentcompare; % get current pair
roi_1_val = handles.pairs_cells1(tt); % cell 1
roi_2_val = handles.pairs_cells2(tt); % cell 2
m = handles.mergecells; % get all the merge pairs in mergelist
m = [m; {mat2str([roi_1_val roi_2_val])} ]; % add current pair
set(handles.mergelist, 'String', m, 'Value', length(m)) % update mergelist
handles.mergecells=m; % update matrix

% go to the next pair (if this was not the last pair)
if handles.pairs_currentcompare ~= length(handles.pairs_cells1)

    % roi pair to plot
    handles.pairs_currentcompare = handles.pairs_currentcompare + 1;

    % ROI(s) to plot
    idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2

    % plot background (Cn or Mn) and ROI contours (idx)
    plot_spatial_components(handles, idx)

    % plot temporal traces of the first pair
    plot_temporal_traces(handles, idx)

    % plot ROI pair info
    plot_roi_info(handles, idx)
end

% store data
guidata(hObject,handles);


function checkmergetxt_Callback(hObject, eventdata, handles)
% hObject    handle to checkmergetxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of checkmergetxt as text
%        str2double(get(hObject,'String')) returns contents of checkmergetxt as a double


% --- Executes on button press in checkmerge.
function checkmerge_Callback(hObject, eventdata, handles)
% hObject    handle to checkmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get all ROIs
tmp_roi = []; % empty vector
for i = 1:length(handles.mergecells) % for each pair
    tmp_roi = [tmp_roi str2num(cell2mat(handles.mergecells(i)))]; % fill with all ROIs (roi1 and roi2's)
end
tbl = tabulate(tmp_roi);    % get frequency table of all values from min to max (include not-occuring values)
val_repeat = tbl(find(tbl(:,2)>1),1); % get repeated ROIs

% feedback user
if isempty(val_repeat)
    set(handles.checkmergetxt, 'String', sprintf(['No repeated ROIs detected -> \n          merge ROIs']));
else
    set(handles.checkmergetxt, 'String', sprintf(['Repeated ROIs detected -> \nperform multi-merge or delete ROIs:\n' num2str(val_repeat')]));
end

% % check if values are repeated
% tbl = tabulate(tmp); % frequency of values
% if ~isempty(find(tbl(:,2)>1)) % find if values (column 2) are repeated
%     repeat_vals = find(tbl(:,2)>1);
%     tmp = []; % create empty tmp var
%     for k=1:length(repeat_vals)
%         % get the value of these cells
%         if k==1
%             tmp = [ num2str( tbl(repeat_vals(k),1)) ];
%         else
%             tmp = [tmp ',' num2str( tbl(repeat_vals(k),1)) ];
%         end
%     end
%     str = sprintf(['    Repeated ROIs detected -> \nperform multi-merge or delete ROIs: \n' num2str(tmp)]);
%     set(handles.checkmergetxt,'String',str); % plot text
% else
%     set(handles.checkmergetxt, 'String', sprintf(['No repeated ROIs detected -> \n          merge ROIs']));
% end

% store data
guidata(hObject,handles);


% --- Executes on button press in multimerge.
function multimerge_Callback(hObject, eventdata, handles)
% hObject    handle to multimerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get all ROIs
tmp_roi = []; % empty vector
for i = 1:length(handles.mergecells) % for each pair
    tmp_roi = [tmp_roi str2num(cell2mat(handles.mergecells(i)))]; % fill with all ROIs (roi1 and roi2's)
end
tbl = tabulate(tmp_roi);    % get frequency table of all values from min to max (include not-occuring values)
val_repeat = tbl(find(tbl(:,2)>1),1); % get repeated ROIs

% examine how to multimerge
tmp_mergelist = []; % empty vector
tmp_mergelist = cell2struct(handles.mergecells, 'ids', length(handles.mergecells)); % get all pairs -> struct array with field 'ids'
idx_pairstommerge=[]; % empty vector to keep track of idx of mergelist pairs to multimerge
for i = 1:length(val_repeat) % per repeated ROI

    % check how often repeated ROI is found throughout the mergelist
    idx_tommerge = []; % empty vector to store idx of mergepairs that include repeated ROI
    for y = 1:size(tmp_mergelist,1) % per ROI
        if ~isempty(find(str2num(tmp_mergelist(y).ids) == val_repeat(i))) % check if ROI is in pair (doesn't matter position 1 or 2)
            idx_tommerge = [idx_tommerge y]; % if so, add idx of mergelist pair (index, not ROI numbers) to idx_tomerge
        end
    end

    % get absolute ROI numbers
    tmp_tommerge = []; % empty vector to keep track of multimerge ROIs
    for j = 1:length(idx_tommerge) % per mergelist pair (idx)
        tmp_tommerge = [tmp_tommerge str2num(tmp_mergelist(idx_tommerge(j)).ids)];
    end
    handles.multimergeidx(i).ids= unique(tmp_tommerge); % per repeated ROI, store absolute values of other ROIs to merge with

    % keep track of all
    idx_pairstommerge = [idx_pairstommerge idx_tommerge]; % store list in new var that includes all idx to keep and the new mergelist idx
end

% loop through multimerge values and identify the multimerge strings that are the same
idx_removemmerge = []; % empty variable for repeated multimerge strings
for k = 1:length(handles.multimergeidx) % per ROI that needs multimerge
    for q = k:length(handles.multimergeidx)
        if q~=k & isequal(handles.multimergeidx(q).ids, handles.multimergeidx(k).ids) % not check same multimerge string
            idx_removemmerge = [idx_removemmerge k]; % keep track of multimerge values that need to go (duplicates)
        end
    end
end
% delete the same multimerge strings
handles.multimergeidx(unique(idx_removemmerge)) = [];

% get full list of normal merge pairs (excl any multimerge pairs)
tmp_mergelist_array=handles.mergecells; % array with all the mergelist pairs
tmp_mergelist_array(unique(idx_pairstommerge))=[]; % remove all multimerge values

% combine all the multimerge strings with the normal merge pairs
if ~isempty(idx_pairstommerge) % make sure we have multimerge
    for k = 1:length(handles.multimergeidx) % per multimerge string
        tmp_mergelist_array = [tmp_mergelist_array; {mat2str(handles.multimergeidx(k).ids)}]; % updated array with normal and multimerge
    end
end
handles.mergecells = tmp_mergelist_array; % store all the normal merge pairs and the multimerge pairs
set(handles.mergelist, 'String', tmp_mergelist_array, 'Value', 1); % update mergelist with normal merge pairs and multimerge values

% give user feedback of which ROI values are still repeated
tmp_roi = []; % empty vector
for i = 1:length(handles.mergecells) % for each pair
    tmp_roi = [tmp_roi str2num(cell2mat(handles.mergecells(i)))]; % fill with all ROIs (roi1 and roi2's)
end
tbl = tabulate(tmp_roi);    % get frequency table of all values from min to max (include not-occuring values)
val_repeat = tbl(find(tbl(:,2)>1),1); % get repeated ROIs

% feedback user
if isempty(val_repeat)
    set(handles.checkmergetxt, 'String', sprintf(['  Multi-merge done -> \nno more repeated ROIs found']));
else
    set(handles.checkmergetxt, 'String', sprintf(['Multi-merge done -> \nmore repeated ROIs found:\n' num2str(val_repeat')]));
end

% delete rois_bad values if we plotted those
if handles.plot_bad == 1 % plot bad ROIs
    handles.rois_bad = [];
end

% store data
guidata(hObject,handles);


% --- Executes on button press in delete4mmerge.
function delete4mmerge_Callback(hObject, eventdata, handles)
% hObject    handle to delete4mmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store that we plot a merge pair
handles.to_plot = 0;

% remove from mergelist
k1 = get(handles.mergelist, 'Value'); % get location
k2 = get(handles.mergelist, 'String'); % get total list
k2(k1,:) = []; % remove from list
handles.mergecells(k1,:) = []; % update mergecells list
set(handles.mergelist, 'String', k2, 'Value', size(k2,1)); % update mergelist

% % add current pair to mergelist
% tt = handles.pairs_currentcompare; % get current pair
% roi_1_val = handles.pairs_cells1(tt); % cell 1
% roi_2_val = handles.pairs_cells2(tt); % cell 2
% m = handles.mergecells; % get all the merge pairs in mergelist
% m = [m; {mat2str([roi_1_val roi_2_val])} ]; % add current pair
%
% tmp = k2(k1,:);

% plot next mergelist ROIs (unless it's the last pair)
if k1 <= length(k2)

    % ROI(s) to plot
    idx = str2num(cell2mat(k2(k1,:)));
    handles.plot_roi = idx;

    % select the next ROI in the dellist
    set(handles.mergelist, 'String', k2, 'Value', k1);

    % plot background (Cn or Mn) and ROI contours (idx)
    plot_spatial_components(handles, idx)

    % plot temporal traces of the first pair
    plot_temporal_traces(handles, idx)

    % plot ROI pair info
    plot_roi_info(handles, idx)

end

% store data
guidata(hObject,handles);


% --- Executes on selection change in mergelist.
function mergelist_Callback(hObject, eventdata, handles)
% hObject    handle to mergelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mergelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mergelist

% store that we plot a mergepair
handles.to_plot = 0;

% ROI(s) to plot
k1 = get(handles.mergelist,'Value');    % location in list
k2 = get(handles.mergelist,'String');   % all pairs in list
tmp = k2(k1);   % actual roi pair
idx = str2num(cell2mat(tmp)); % values ROI1 and 2

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% plot ROI pair info
plot_roi_info(handles, idx)

% store data
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS


function a_min_Callback(hObject, eventdata, handles)
% hObject    handle to a_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_min as text
%        str2double(get(hObject,'String')) returns contents of a_min as a double


function a_max_Callback(hObject, eventdata, handles)
% hObject    handle to a_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_max as text
%        str2double(get(hObject,'String')) returns contents of a_max as a double


function snr_thres_Callback(hObject, eventdata, handles)
% hObject    handle to snr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr_thres as text
%        str2double(get(hObject,'String')) returns contents of snr_thres as a double


function dist_thres_Callback(hObject, eventdata, handles)
% hObject    handle to dist_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist_thres as text
%        str2double(get(hObject,'String')) returns contents of dist_thres as a double
% handles.distthr =  str2double(get(handles.dist_thres,'String'))


function corr_thres_Callback(hObject, eventdata, handles)
% hObject    handle to corr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corr_thres as text
%        str2double(get(hObject,'String')) returns contents of corr_thres as a double
% handles.corrthr =  str2double(get(handles.corr_thres,'String'))


% --- Executes on button press in check_plot_bad.
function check_plot_bad_Callback(hObject, eventdata, handles)
% hObject    handle to check_plot_bad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_plot_bad

% check if we allow updating plotting (i.e., have not loaded data yet)
if handles.plot_bad_update == 0 % prevent to update
    set(handles.check_plot_bad, 'Value', handles.plot_bad); % 0=good ROIs, 1=bad ROIs

    % update user
    set(handles.user_alert, 'String', sprintf('You already loaded data, cannot plot bad ROIs!'));
else

    % update user
    set(handles.user_alert, 'String', sprintf('Important: plot bad ROIs of raw output only once!'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING


function contour_thres_Callback(hObject, eventdata, handles)
% hObject    handle to contour_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contour_thres as text
%        str2double(get(hObject,'String')) returns contents of contour_thres as a double

% ROI(s) to plot
if handles.to_plot == 0; % plot mergepair
    idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2
elseif handles.to_plot == 1; % plot single ROI
    idx = handles.plot_roi; % plot roi
elseif handles.to_plot == 2; % plot manual ROI pair
    idx = handles.manualpair;
end

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% store data
guidata(hObject,handles);


function scale_cn_Callback(hObject, eventdata, handles)
% hObject    handle to scale_cn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale_cn as text
%        str2double(get(hObject,'String')) returns contents of scale_cn as a double

% ROI(s) to plot
if handles.to_plot == 0; % plot mergepair
    idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2
elseif handles.to_plot == 1; % plot single ROI
    idx = handles.plot_roi; % plot roi
elseif handles.to_plot == 2; % plot manual ROI pair
    idx = handles.manualpair;
end

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% % plot temporal traces of the first pair
% plot_temporal_traces(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on selection change in cn_mn_plot.
function cn_mn_plot_Callback(hObject, eventdata, handles)
% hObject    handle to cn_mn_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cn_mn_plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cn_mn_plot

% ROI(s) to plot
if handles.to_plot == 0; % plot mergepair
    idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2
elseif handles.to_plot == 1; % plot single ROI
    idx = handles.plot_roi; % plot roi
elseif handles.to_plot == 2; % plot manual ROI pair
    idx = handles.manualpair;
end

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% % plot temporal traces of the first pair
% plot_temporal_traces(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on selection change in craw_c.
function craw_c_Callback(hObject, eventdata, handles)
% hObject    handle to craw_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns craw_c contents as cell array
%        contents{get(hObject,'Value')} returns selected item from craw_c

% ROI(s) to plot
if handles.to_plot == 0; % plot mergepair
    idx = [handles.pairs_cells1(handles.pairs_currentcompare) handles.pairs_cells2(handles.pairs_currentcompare)]; % ROI 1 and ROI 2
elseif handles.to_plot == 1; % plot single ROI
    idx = handles.plot_roi; % plot roi
elseif handles.to_plot == 2; % plot manual ROI pair
    idx = handles.manualpair;
end

% % plot background (Cn or Mn) and ROI contours (idx)
% plot_spatial_components(handles, idx)

% plot temporal traces of the first pair
plot_temporal_traces(handles, idx)

% store data
guidata(hObject,handles);


% --- Executes on button press in plot_rois.
function plot_rois_Callback(hObject, eventdata, handles)
% hObject    handle to plot_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_rois

% update user
set(handles.user_alert, 'String', sprintf('Plotting all ROIs, please wait'));

% ROI(s) to plot
idx = 1:size(handles.updatedCraw,1); % all ROIs

% plot background (Cn or Mn) and ROI contours (idx)
plot_spatial_components(handles, idx)

% turn off radio button
set(handles.plot_rois, 'Value', 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SINGLE ROI DATA


function roi_idx_Callback(hObject, eventdata, handles)
% hObject    handle to roi_idx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi_idx as text
%        str2double(get(hObject,'String')) returns contents of roi_idx as a double


function n_rois_Callback(hObject, eventdata, handles)
% hObject    handle to n_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_rois as text
%        str2double(get(hObject,'String')) returns contents of n_rois as a double


function roi_snr_Callback(hObject, eventdata, handles)
% hObject    handle to roi_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi_snr as text
%        str2double(get(hObject,'String')) returns contents of roi_snr as a double


function roi_size_Callback(hObject, eventdata, handles)
% hObject    handle to roi_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi_size as text
%        str2double(get(hObject,'String')) returns contents of roi_size as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROI PAIR DATA


function cell_pair1_Callback(hObject, eventdata, handles)
% hObject    handle to cell_pair1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_pair1 as text
%        str2double(get(hObject,'String')) returns contents of cell_pair1 as a double


function cell_pair_dist_Callback(hObject, eventdata, handles)
% hObject    handle to cell_pair_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_pair_dist as text
%        str2double(get(hObject,'String')) returns contents of cell_pair_dist as a double


function cell_pair2_Callback(hObject, eventdata, handles)
% hObject    handle to cell_pair2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_pair2 as text
%        str2double(get(hObject,'String')) returns contents of cell_pair2 as a double


function cell_pair_corr_Callback(hObject, eventdata, handles)
% hObject    handle to cell_pair_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_pair_corr as text
%        str2double(get(hObject,'String')) returns contents of cell_pair_corr as a double


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE PLOT PANELS


% --- Executes during object creation, after setting all properties.
function distcorr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distcorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate distcorr


% --- Executes during object creation, after setting all properties.
function histsnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histsnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate histsnr


% --- Executes during object creation, after setting all properties.
function ax1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax1


% --- Executes during object creation, after setting all properties.
function ax2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE ALL OTHER OBJECTS


% --- Executes during object creation, after setting all properties.
function roi1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function roi2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function dist_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function corr_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function contour_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contour_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function snr_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function deletelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deletelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function celllist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to celllist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function add2dellist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to add2dellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function slidertoaddtodellist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidertoaddtodellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function mergelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mergelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');mlist
end


% --- Executes during object creation, after setting all properties.
function checkmerge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function checkmergetxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkmergetxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function id_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function scale_cn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale_cn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cn_mn_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cn_mn_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cell_pair1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_pair1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cell_pair_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_pair_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cell_pair2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_pair2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cell_pair_corr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_pair_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function roi_idx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_idx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function roi_snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on mergelist and none of its controls.
function mergelist_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to mergelist (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function craw_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to craw_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function n_rois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function user_alert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_alert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function a_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function a_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function roi_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function n_pair_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_sizehisto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_sizehisto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function histsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate histsize


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL FUNCTIONS


function [snr] = get_SNR(dataC, dataCraw)
% estimate SNR based on peaks / std(noise)
% find location peaks based on C, value based on Craw
% noise is all Craw data without peaks
% inputs:
%   dataC = denoised trace to calculate peaks, T*1 vector, calcium trace
%   dataCraw = raw trace to calculate noise, T*1 vector, calcium trace
%   perc = percentile to use to estimate baseline

% fixed vars
windowCraw = 5; % window before peak in C to find peak in Craw data
windowCrawRemove = 2; % window around Craw peak to remove for baseline
baselinePercentile = 90; % percentile baseline to infer noise

% find peaks and baseline to calculate SNR
for k = 1:size(dataC,1)    % per roi

    % get baseline data Craw
    baseline = dataCraw(k,:);

    % peaks in C
    [~, loc] = findpeaks(dataC(k,:));

    % get window to find peaks in Craw
    locOn = loc - windowCraw; % go back in time to find max peak in Craw
    locOff = loc;

    % make sure locOn is above frame 1
    for kk = 1:length(locOn) % per peak 
        if locOn(kk) < 1
            locOn(kk) = 1;
        end
    end

    % find peak values in Craw and remove from baseline
    for kk = 1:length(locOn) % per peak
        [max_val(kk), max_idx(kk)] = max(dataCraw(k, locOn(kk):locOff(kk))); % get max value and its location
        abs_max_idx(kk) = max_idx(kk)+locOn(kk)-1; % find absolute location
        if abs_max_idx(kk)-windowCrawRemove < 1 % value might be frame 1
            baseline(1, 1:abs_max_idx(kk)+2) = nan(1,length(1:abs_max_idx(kk)+2));
        else
            baseline(1, abs_max_idx(kk)-windowCrawRemove : abs_max_idx(kk)+windowCrawRemove) = nan(1,5);
        end
    end
    
    % calculate SNR
    if exist('max_val', 'var') % check if we have at least 1 max_val
%         snr(k) = nanmedian(max_val) / nanstd(baseline);
        snr(k) = nanmedian(max_val) / prctile(baseline, baselinePercentile);
    else
        snr(k) = nan;
    end

    % clear vars
    clear loc locOn locOff baseline max_val max_idx abs_max_idx
end

% % find peaks and baseline to calculate SNR
% for k = 1:size(dataC,1)    % per roi
%     % get baseline data
%     baseline = dataCraw(k,:);
%     % peaks in C
%     [~, loc] = findpeaks(dataC(k,:)); % , 3
%     % find peak values in Craw
%     pk = dataCraw(k,loc);
%     % baseline in Craw
%     baseline(loc) = [];
%     % SNR
%     snr(k) = median(pk) / std(baseline);
% 
%     clear pk loc baseline
% end


function plot_temporal_traces(handles, idx)
% plot temporal components in ax1

% fixed vars
fraction_yrange = .05;   % fraction of data to plot above/below figure

% check if we want to plot Craw or C
plot_ccraw = get(handles.craw_c, 'Value');

% plot Craw or C
min_val = []; % keep track of min val y axis
max_val = []; % keep track of max val y axis
axes(handles.ax1);
if plot_ccraw == 1 % C_raw

    % identify 1 ROI, 2 ROIs, or anything else (multi-merge)
    if length(idx) == 1 % one ROI
        plot(handles.ax1, handles.updatedCraw(idx, :), 'Color', handles.colors(1,:));
        min_val = min(handles.updatedCraw(idx, :));
        max_val = max(handles.updatedCraw(idx, :));
    elseif length(idx) == 2 % 2 ROIs
        % plot per ROI
        for k = 1:length(idx)
            plot(handles.ax1, handles.updatedCraw(idx(k), :), 'Color', handles.colors(k+1,:));
            hold(handles.ax1,'on')
            min_val = [min_val min(handles.updatedCraw(idx(k), :))];
            max_val = [max_val max(handles.updatedCraw(idx(k), :))];
        end
        hold(handles.ax1,'off')
    else % many ROIs multi-merge
        % plot per ROI
        cmap = parula(length(idx));
        for k = 1:length(idx)
            plot(handles.ax1, handles.updatedCraw(idx(k), :), 'Color', cmap(k,:));
            hold(handles.ax1,'on')
            min_val = [min_val min(handles.updatedCraw(idx(k), :))];
            max_val = [max_val max(handles.updatedCraw(idx(k), :))];
        end
        hold(handles.ax1,'off')
    end

    % define y lims
    try
        change_y = (max(max_val) - min(min_val)) * fraction_yrange;
        ylim([min(min_val)-change_y max(max_val)+change_y])
    catch end % in case of a rois_bad, signal can be 0

    % define x lims
    xlim([1 size(handles.updatedCraw,2)])

else % plot C

    % identify 1 ROI, 2 ROIs, or anything else (multi-merge)
    if length(idx) == 1 % one ROI
        plot(handles.ax1, handles.updatedC(idx, :), 'Color', handles.colors(1,:));
        min_val = min(handles.updatedC(idx, :));
        max_val = max(handles.updatedC(idx, :));
    elseif length(idx) == 2 % 2 ROIs
        % plot per ROI
        for k = 1:length(idx)
            plot(handles.ax1, handles.updatedC(idx(k), :), 'Color', handles.colors(k+1,:));
            hold(handles.ax1,'on')
            min_val = [min_val min(handles.updatedC(idx(k), :))];
            max_val = [max_val max(handles.updatedC(idx(k), :))];
        end
        hold(handles.ax1,'off')
    else % many ROIs multi-merge
        % plot per ROI
        cmap = parula(length(idx));
        for k = 1:length(idx)
            plot(handles.ax1, handles.updatedC(idx(k), :), 'Color', cmap(k,:));
            hold(handles.ax1,'on')
            min_val = [min_val min(handles.updatedC(idx(k), :))];
            max_val = [max_val max(handles.updatedC(idx(k), :))];
        end
        hold(handles.ax1,'off')
    end

    % define y lims
    try
        change_y = (max(max_val) - min(min_val)) * fraction_yrange;
        ylim([min(min_val)-change_y max(max_val)+change_y])
    catch end

    % define x lims
    xlim([1 size(handles.updatedCraw,2)])

end


function plot_spatial_components(handles, idx)
% plot background (Cn or Mn) and spatial components in ax2

%///////////////////////////////////// PLOT BACKGROUND
% get figure
axes(handles.ax2);

% get c axis scale
c_lim_reduction =  str2num(get(handles.scale_cn,'String'));

% get what to plot (Cn or Mn) and plot
plot_cn_mn = get(handles.cn_mn_plot,'Value');

% plot figure
if plot_cn_mn == 1 % plot Cn
    imagesc(fliplr(rot90(handles.rawdata1.results.Cn, -1)))
    tmp_max = ( min(min(handles.rawdata1.results.Cn)) + ((max(max(handles.rawdata1.results.Cn)) - min(min(handles.rawdata1.results.Cn))) / c_lim_reduction) );
    clim([ min(min(handles.rawdata1.results.Cn)) + ((tmp_max - min(min(handles.rawdata1.results.Cn))) / 2) tmp_max]);
elseif plot_cn_mn == 2 % plot Mn
    imagesc(fliplr(rot90(handles.rawdata1.results.Mn, -1)))
    clim([min(min(handles.rawdata1.results.Mn)) ...
        ( min(min(handles.rawdata1.results.Mn)) + ((max(max(handles.rawdata1.results.Mn)) - min(min(handles.rawdata1.results.Mn))) / c_lim_reduction) )]);
end

% change colormap
colormap gray;

%///////////////////////////////////// PLOT ROIs
% get figure
hold(handles.ax2,'on')

% plot per ROI
for k = 1:length(idx)
    % extract the spatial component and reshape it
%     roiSpatial = full(handles.updatedA(:, idx(k))); % convert sparse column to full
%     roiImage = reshape(roiSpatial, handles.rawdata1.results.options.d1, handles.rawdata1.results.options.d2); % reshape to 2D
    roiImage = squeeze(handles.Afull(idx(k),:,:));

    % plot contour
    % [~, h] = contour(roiImage, [str2num(get(handles.contour_thres,'String')) str2num(get(handles.contour_thres,'String'))], 'LineWidth', 1.5); % adjust threshold as needed
    tmp_thres = min(min(roiImage)) + ((max(max(roiImage)) - min(min(roiImage))) * (1 - str2num(get(handles.contour_thres,'String'))) ); % contour_thres is fraction of entire ROI

    % identify 1 ROI, 2 ROIs, or anything else (multi-merge)
    if length(idx) == 1 % one ROI = del list
        contour(roiImage, [tmp_thres tmp_thres], 'LineWidth', 1.5, 'LineColor', handles.colors(1,:)); % adjust threshold as needed
    elseif length(idx) == 2 % two ROIs = merge list
        contour(roiImage, [tmp_thres tmp_thres], 'LineWidth', 1.5, 'LineColor', handles.colors(k+1,:)); % adjust threshold as needed
    else % multi-merge
        cmap = parula(length(idx));
        contour(roiImage, [tmp_thres tmp_thres], 'LineWidth', 1.5, 'LineColor', cmap(k,:)); % adjust threshold as needed
    end
end

% remove labels
set(handles.ax2, 'XTick', [], 'YTick', []);
hold(handles.ax2,'off')


function plot_roi_info(handles, idx)
% plot info regarding ROIs
% 1 ROI: update roi index and SNR
% 2 ROIs: roi indices, dist and corr, and dist vs corr plot (ax3)

% figure out how many ROIs
if length(idx) == 1 % 1 roi
    % display roi index and SNR value
    set(handles.roi_idx, 'String', sprintf('%d', idx));
    set(handles.roi_snr, 'String', sprintf('%.1f', handles.SNR(idx)));
    set(handles.roi_size, 'String', sprintf('%.0f', handles.pixels_A(idx)));

else % more ROIs
    % plot cell identity
    set(handles.cell_pair1, 'String', sprintf('%d', idx(1)));
    set(handles.cell_pair2, 'String', sprintf('%d', idx(2)));

    % if we plot a manual pair, sort idx to make high value roi 1
    if handles.to_plot == 2 % manual merge pair
        idx = sort(idx, 'descend');
    end

    % plot cell pair info
    set(handles.cell_pair_dist, 'String', sprintf('%.1f', handles.dist_mat(idx(1), idx(2))));
    set(handles.cell_pair_corr, 'String', sprintf('%.1f', handles.corr_mat(idx(1), idx(2))));

    % plot distance vs crosscorr: all data
    plot(handles.distcorr, handles.dist, handles.corr, '.', 'Color', handles.color_plots(1,:), 'LineStyle', 'none', 'MarkerSize', 10);
    hold(handles.distcorr,'on')

    % plot distance vs crosscorr: current pair
    plot(handles.distcorr, handles.dist_mat(idx(1), idx(2)), handles.corr_mat(idx(1), idx(2)), '.', 'Color', handles.color_plots(2,:), 'LineStyle', 'none', 'MarkerSize', 20)
    hold(handles.distcorr,'off')
end


function output = load_cnmf_hdf5(filename1, filepath1, handles);
% load raw hdf5 output of caiman (cnmf)
% store data so that the GUI can use it
% keep a .raw with all the data (good/bad rois, settings, etc)

% identify structure and database
metadata = h5info([filepath1 filename1]); % read file
data_nm = metadata.Groups(1).Name;
options_nm = metadata.Groups(2).Name;

% create struct to load data
output = struct;

% load the data to generate sparse matrix for contours (A)
data = h5read([filepath1 filename1], '/estimates/A/data');
indices = h5read([filepath1 filename1], '/estimates/A/indices') + 1; % MATLAB is 1-based
indptr = h5read([filepath1 filename1], '/estimates/A/indptr') + 1; % MATLAB is 1-based
shape = h5read([filepath1 filename1], '/estimates/A/shape');

% reconstruct the sparse matrix
nRows = double(shape(1));
nCols = double(shape(2));
row_indices = [];
col_indices = [];
values = [];

% generate values for matrix
for col = 1:nCols
    start_idx = indptr(col);
    end_idx = indptr(col+1) - 1;
    row_indices = [row_indices; indices(start_idx:end_idx)];
    col_indices = [col_indices; double(col) * ones(end_idx - start_idx + 1, 1)]; % explicitly cast col to double
    values = [values; data(start_idx:end_idx)];
end

% ensure values are double as well
row_indices = double(row_indices);
col_indices = double(col_indices);
values = double(values);

% create the sparse matrix
A = sparse(row_indices, col_indices, values, nRows, nCols);

% find caiman data
% C, Cn, F_dff, R, S, SNR_comp, YrA, bl, c1, g, idx_components,
% idx_components_bad, lam, neurons_sn, r_values, data_nm
C_raw = h5read([filepath1 filename1], [data_nm '/F_dff'])';
C = h5read([filepath1 filename1], [data_nm '/C'])';
S = h5read([filepath1 filename1], [data_nm '/S'])';
spat_bgnd_comp = h5read([filepath1 filename1], [data_nm '/b']);
temp_bgnd_comp = h5read([filepath1 filename1], [data_nm '/f']);
rois_good = h5read([filepath1 filename1], [data_nm '/idx_components']) +1; % MATLAB is 1-based
rois_bad = h5read([filepath1 filename1], [data_nm '/idx_components_bad']) +1; % MATLAB is 1-based
rois_bad = rois_bad'; %
try % not everyone has these vars
    residuals = h5read([filepath1 filename1], [data_nm '/YrA'])';
    tmp = h5read([filepath1 filename1], [data_nm '/neurons_sn']); % get SNR
    for k=1:length(tmp)
        SNR(k,1) = str2num(tmp(k));
    end
catch end

% store raw rois
output.raw.C_raw = C_raw;
output.raw.C = C;
output.raw.S = S;
output.raw.A = A;
output.raw.spat_bgnd_comp = spat_bgnd_comp;
output.raw.temp_bgnd_comp = temp_bgnd_comp;
try % not everyone has these vars
    output.raw.residuals = residuals;
    output.raw.SNR_cnmf = SNR;
catch end
output.raw.rois_good = rois_good;
output.raw.rois_bad = rois_bad;

% try to include motion correction shifts
try % not everyone has the (pw)rigid movements
    output.raw_correlations = h5read([filepath1 filename1], [options_nm '/raw_corrs']);
    output.rigid_correlations = h5read([filepath1 filename1], [options_nm '/rigid_corrs']);
    output.rigid_shifts = h5read([filepath1 filename1], [options_nm '/rigid_shifts']);
    output.pwrigid_shifts = h5read([filepath1 filename1], [options_nm '/pwrigid_shifts']);
catch end

% store good or all ROIs -> if all, the bad ROIs go into delete list
if handles.plot_bad == 0 % good ROIs
    % store only the good rois
    output.C_raw = C_raw(rois_good, :);
    output.C = C(rois_good, :);
    output.S = S(rois_good, :);
    output.A = A(:, rois_good);
    try % not everyone has these vars
        output.residuals = residuals(rois_good, :);
        output.SNR_cnmf = SNR(rois_good, :);
    catch end
elseif handles.plot_bad == 1 % all ROIs (good and bad)
    % store all the ROIs
    output.C_raw = C_raw(:, :);
    output.C = C(:, :);
    output.S = S(:, :);
    output.A = A(:, :);
    try % not everyone has these vars
        output.residuals = residuals(:, :);
        output.SNR_cnmf = SNR(:, :);
    catch end
    output.rois_bad = rois_bad;
end

% store images
output.Cn = h5read([filepath1 filename1], [data_nm '/Cn']);   % correlation image
try % not everyone has Mn
    output.Mn = h5read([filepath1 filename1], [data_nm '/Mn']);   % max projection
    output.Cn_ori = h5read([filepath1 filename1], [data_nm '/Cn_ori']);   % correlation image original
    output.Mn_ori = h5read([filepath1 filename1], [data_nm '/Mn_ori']);   % max projection original
catch end

% store caiman options
output.options.fps = h5read([filepath1 filename1], [options_nm '/data/fr']);
output.options.dims = h5read([filepath1 filename1], [options_nm '/data/dims']);
output.options.d1 = output.options.dims(1);
output.options.d2 = output.options.dims(2);
output.options.decay_time = h5read([filepath1 filename1], [options_nm '/data/decay_time']);
output.options.fnames = h5read([filepath1 filename1], [options_nm '/data/fnames']);
output.options.dxy = h5read([filepath1 filename1], [options_nm '/data/dxy']);


function [center, pixels, A_full] = get_shape_size_A(handles)
% get shape of all ROIs, find center of mass, and find number of pixels

% get shape of ROIs A
A_full = reshape(full(handles.updatedA)', size(handles.updatedCraw,1), handles.rawdata1.results.options.d1, handles.rawdata1.results.options.d2);

% loop through cells to calculate center of mass and find size A
center = zeros(size(handles.updatedCraw,1),2);
for i = 1:size(handles.updatedCraw,1)
    % COM
    [tmp_r,tmp_c] = find(squeeze(A_full(i,:,:))); % get all rows and columns
    center(i,:) = [mean(tmp_r) mean(tmp_c)]; % find mean of row and column
    % size A
    pixels(i) = length(find(A_full(i,:,:) ~= 0)); % get number of pixels of A
end


function plot_histo_snr(handles)
% plot histogram of SNR

% histsnr
axes(handles.histsnr);
histogram(handles.SNR, 'FaceColor', handles.color_plots(1,:))


function plot_histo_size(handles)
% plot histogram of size

%histsize
axes(handles.histsize);
h = histogram(handles.pixels_A, 'FaceColor', handles.color_plots(1,:));
xlim([0 (h.BinEdges(end) * 1.15)]) % change x lims
hold(handles.histsize,'on')
line([str2num(get(handles.a_min,'String')) str2num(get(handles.a_min,'String'))], [0 max(h.Values)], 'LineWidth', 1, 'Color', handles.color_plots(2,:))% plot manual A_min
line([str2num(get(handles.a_max,'String')) str2num(get(handles.a_max,'String'))], [0 max(h.Values)], 'LineWidth', 1, 'Color', handles.color_plots(2,:))% plot manual A_max

hold(handles.histsize,'off')
