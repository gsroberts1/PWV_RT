function varargout = PWV_RT(varargin)
    % PWV_RT MATLAB code for PWV_RT.fig
    %      PWV_RT, by itself, creates a new PWV_RT or raises the existing
    %      singleton*.
    %
    %      H = PWV_RT returns the handle to a new PWV_RT or the handle to
    %      the existing singleton*.
    %
    %      PWV_RT('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in PWV_RT.M with the given input arguments.
    %
    %      PWV_RT('Property','Value',...) creates a new PWV_RT or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before PWV_RT_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to PWV_RT_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help PWV_RT

    % Last Modified by GUIDE v2.5 20-Jun-2021 22:27:30
    % Developed by Grant S Roberts, University of Wisconsin-Madison, 2019
    
    
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @PWV_RT_OpeningFcn, ...
                       'gui_OutputFcn',  @PWV_RT_OutputFcn, ...
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


function PWV_RT_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to PWV_RT (see VARARGIN)

    % Choose default command line output for PWV_RT
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % Get anatomical, pc, and bSSFP data from LoadPWV GUI
    handles.pcDatasets(1).ROI = []; %initialize ROI field
    handles.pcDatasets(1).ROIdata = []; %flow smoothed with Gauss.
    
    handles.global.startAnalyzing = 0; %flag to begin PWV calculations
    handles.global.totalROIs = 0;
    handles.global.pcIter = 1;

    set(handles.load2DPCbutton,'Enable','off');
    set(handles.pcPlanePopup,'Enable','off');
    set(handles.pcDatasetPopup,'Enable','off');
    set(handles.drawROIbutton,'Enable','off');
    set(handles.loadROIbutton,'Enable','off');
    set(handles.pcSlider,'Enable','off');
    set(handles.ttpointRadio,'Value',1); 
    set(handles.ttpointRadio,'Enable','off');
    set(handles.ttuRadio,'Value',1);
    set(handles.ttuRadio,'Enable','off');
    set(handles.ttfRadio,'Value',1);
    set(handles.ttfRadio,'Enable','off');
    set(handles.xcorrRadio,'Value',1);
    set(handles.xcorrRadio,'Enable','off');
    set(handles.exportAnalysisButton,'Enable','off');
    
    guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = PWV_RT_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;


    
%%%%%%%%%%%% LOAD 2DPC PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PLANE PLOT - CREATE FUNCTION
function pcPlanePlot_CreateFcn(hObject, eventdata, handles)


% --- LOAD CENTERLINE DATA - CALLBACK
function loadCLpush_Callback(hObject, eventdata, handles)
    [clFile, clDir] = uigetfile({'*.mat;','Useable Files (*.mat)';
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select the centerline dataset (anatCLdataset.mat)');
    load([clDir clFile]);
    handles.centerline = anatCLdataset;
    
    cd(clDir);
    handles.global.homeDir = clDir;
    set(handles.load2DPCbutton,'Enable','on');
    guidata(hObject, handles);
    
    
% --- LOAD 2DPC DATASETS - CALLBACK
function load2DPCbutton_Callback(hObject, eventdata, handles)
    [pcFile, pcDir] = uigetfile({'*.dcm;*.dat;*.mat','Useable Files (*.dcm,*.dat,*.mat)';
       '*.dcm',  'DICOM files (*.dcm)'; ...
       '*.dat',  'DAT-files (*.dat)'; ...
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE 2DPC file in the dataset');
    pcIter = handles.global.pcIter;
    [~,~,extension] = fileparts(pcFile);
    dirInfo = dir(fullfile(pcDir,['*' extension]));
    if isequal(extension,'.dcm') %if our extension is a dicom file
        handles.pcDatasets(pcIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name)); %get dicom metadata (from 1st dicom)
        for i=1:length(dirInfo)
            hold(:,:,i) = single(dicomread(fullfile(pcDir,dirInfo(i).name))); %read dicoms and cast to single
        end  
        mag = hold(:,:,floor(length(dirInfo)/2)+1:end); %magnitude is last half
        v = hold(:,:,1:floor(length(dirInfo)/2)); %velocity is first half of images
        MAG = mean(mag,3); %time-averaged magnitude
        VMEAN = mean(v,3); %time-averaged velocity
        CD = MAG.*sin( pi/2*abs(VMEAN)/max(VMEAN(:)) );
        
        handles.pcDatasets(pcIter).Images.MAG = MAG;
        handles.pcDatasets(pcIter).Images.CD = CD;
        handles.pcDatasets(pcIter).Images.V = VMEAN;
        handles.pcDatasets(pcIter).Images.mag = mag;
        handles.pcDatasets(pcIter).Images.v = v; 
        handles.pcDatasets(pcIter).Names = ['Plane' num2str(pcIter)];
    elseif isequal(extension,'.dat')
        fid = fopen([pcDir filesep 'pcvipr_header.txt'], 'r'); %open header
        dataArray = textscan(fid,'%s%s%[^\n\r]','Delimiter',' ', ...
            'MultipleDelimsAsOne',true,'ReturnOnError',false); %parse header info
        fclose(fid);
        dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:),'UniformOutput',false);
        pcviprHeader = cell2struct(dataArray{1,2}(:),dataArray{1,1}(:),1); %turn to structure
        handles.pcDatasets(pcIter).Info = pcviprHeader; %add pcvipr header to handles
        resx = pcviprHeader.matrixx; %resolution in x
        resy = pcviprHeader.matrixy; %resolution in y
        nframes = pcviprHeader.frames; %number of cardiac frames
        MAG = load_dat(fullfile(pcDir,'MAG.dat'),[resx resy]); %Average magnitude
        CD = load_dat(fullfile(pcDir,'CD.dat'),[resx resy]); %Average complex difference
        VMEAN = load_dat(fullfile(pcDir,'comp_vd_3.dat'),[resx resy]); %Average velocity

        % Initialize data time-resolved data arrays
        v = zeros(resx,resy,nframes); %Time-resolved velocity 
        for j = 1:nframes  %velocity is placed in v3 for 2D (through-plane)
            v(:,:,j) = load_dat(fullfile(pcDir,[filesep 'ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
        end
        
        handles.pcDatasets(pcIter).Images.MAG = flipud(MAG);
        handles.pcDatasets(pcIter).Images.CD = flipud(CD);
        handles.pcDatasets(pcIter).Images.V = flipud(VMEAN);
        handles.pcDatasets(pcIter).Images.v = flipud(v); 
        handles.pcDatasets(pcIter).Names = ['Plane' num2str(pcIter)];
    else %if a single matlab file (with all images)
        hold = load([pcDir pcFile]);
        planeName = fieldnames(hold);
        planeName = planeName{1};
        handles.pcDatasets(pcIter).Info = hold.(planeName).Info;
        images = hold.(planeName).Images;
        images = single(images);
        
        handles.pcDatasets(pcIter).Images.MAG = mean(images(:,:,:,1),3);
        handles.pcDatasets(pcIter).Images.CD = mean(images(:,:,:,2),3);
        handles.pcDatasets(pcIter).Images.V = mean(images(:,:,:,3),3);
        handles.pcDatasets(pcIter).Images.mag = images(:,:,:,1);
        handles.pcDatasets(pcIter).Images.cd = images(:,:,:,2);
        handles.pcDatasets(pcIter).Images.v = images(:,:,:,3);
        handles.pcDatasets(pcIter).Names = planeName;
    end
    
    handles.global.pcIter = handles.global.pcIter + 1;
    
    set(handles.pcPlanePopup,'Enable','on');
    set(handles.pcDatasetPopup,'Enable','on');
    set(handles.drawROIbutton,'Enable','on');
    set(handles.loadROIbutton,'Enable','on');
    set(handles.pcSlider,'Enable','on');
    set(handles.pcPlanePopup,'String',{handles.pcDatasets.Names}); %list of all planes (AAo, AbdAo, etc.)
    set(handles.pcDatasetPopup,'String',fieldnames(handles.pcDatasets(pcIter).Images)); %list of all datasets (CD, MAG, v, etc.)
    guidata(hObject, handles);
    updatePCImages(handles);
    
    
% --- PLANE DROPDOWN - CALLBACK
function pcPlanePopup_Callback(hObject, eventdata, handles)   
    updatePCImages(handles); %update images on PC plot anytime we click on a new plane

% --- PLANE DROPDOWN - CREATE FUNCTION
function pcPlanePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- DATASET DROPDOWN - CALLBACK
function pcDatasetPopup_Callback(hObject, eventdata, handles)
    updatePCImages(handles); %update images on PC plot anytime we click on a new dataset

% --- DATASET DROPDOWN - CREATE FUNCTION
function pcDatasetPopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- DRAWROI BUTTON - CALLBACK
function drawROIbutton_Callback(hObject, eventdata, handles)
    axes(handles.pcPlanePlot); %make sure we do this on top PC plot 

    set(handles.loadCLpush,'Enable','off');
    set(handles.load2DPCbutton,'Enable','off');
    set(handles.pcPlanePopup,'Enable','off'); %make it so we can't select another plane
    set(handles.pcDatasetPopup,'Enable','off'); %make it so we can't select another plane
    set(handles.drawROIbutton,'Enable','off'); %make it so we can't draw a second ROI
    
    planeNum = get(handles.pcPlanePopup,'Value'); %get current plane of interest (eg AAo)
    mydlg = warndlg('Press enter when the ROI is set'); %open dialog warning 
    waitfor(mydlg); %MUST PRESS ENTER TO PROCEED
    circle = drawcircle('FaceAlpha',0.1,'Color','g','LineWidth',1,'Deletable',0); %draw circle on PC image
    while true
        w = waitforbuttonpress; %wait for enter push ...
        switch w 
            case 1 % if it was a keyboard press.
            key = get(gcf,'currentcharacter'); %get key that was pressed
                switch key
                    case 27 % escape key
                        set(handles.loadCLpush,'Enable','on');
                        set(handles.load2DPCbutton,'Enable','on');
                        set(handles.pcPlanePopup,'Enable','on'); %make it so we can't select another plane
                        set(handles.pcDatasetPopup,'Enable','on'); %make it so we can't select another plane
                        set(handles.drawROIbutton,'Enable','on');
                        broke = 1;
                        break % break out of the while loop
                    case 13 % 13 is the enter/return key 
                        circle.InteractionsAllowed = 'none'; %freeze circle
                        broke = 0;
                        break
                    otherwise 
                        %wait for a different command
                end
       end
    end
    
    if ~broke
        handles.pcDatasets(planeNum).ROI = circle; %temporarily hold circle (deleted once gone)
    else 
        handles.pcDatasets(planeNum).ROI = [];
    end 

    guidata(hObject,handles);
    updatePCImages(handles); 

    
% --- LOAD ROI BUTTON - CALLBACK
function loadROIbutton_Callback(hObject, eventdata, handles)
    planeNum = get(handles.pcPlanePopup,'Value'); %get current plane (eg AAo)
    handles.global.totalROIs = handles.global.totalROIs + 1; %add 1 to total ROI count

    v = handles.pcDatasets(planeNum).Images.v; %grab time-resolved velocity
    circle = handles.pcDatasets(planeNum).ROI; %pull circle data from handles
    radius = circle.Radius; %get radius of circle
    center = round(circle.Center); %get center coordinates

    [X,Y] = ndgrid(1:size(v,1),1:size(v,2));
    X = X-center(2); %shift coordinate grid
    Y = Y-center(1);
    roiMask = sqrt(X.^2 + Y.^2)<=radius; %anything outside radius is ignored

    %%% Create Linear Interpolated Data
    if isfield(handles.pcDatasets(planeNum).Info,'matrixx') %if radial data (pcvipr recon)
        matrixx = handles.pcDatasets(planeNum).Info.matrixx; %matrix size in x dimension
        fovx = handles.pcDatasets(planeNum).Info.fovx;  %field of view (mm)
        xres = fovx/matrixx; %resolution (mm). ASSUMED TO BE SAME IN Y DIMENSION
        frames = handles.pcDatasets(planeNum).Info.frames;
        timeres = handles.pcDatasets(planeNum).Info.timeres; %temporal resolution (ms)
%         frames = handles.pcDatasets(1).Info.frames;
%         timeres = handles.pcDatasets(1).Info.timeres;
    else 
        xres = handles.pcDatasets(planeNum).Info.PixelSpacing(1); %resolution (mm) ASSUMED SAME IN Y DIM
        rrInterval = handles.pcDatasets(planeNum).Info.NominalInterval; %average RR interval (ms)
%         rrInterval = handles.pcDatasets(1).Info.NominalInterval; %avg RR int. (ms) TAKE FIRST TO MATCH TEMP RES
        frames = handles.pcDatasets(planeNum).Info.CardiacNumberOfImages; %number of cardiac frames
        timeres = rrInterval/frames; %temporal resolution (ms)
    end 

    area = sum(roiMask(:))*(xres)^2; %ROI area (mm^2)
    for i=1:frames
        vTemp = v(:,:,i); %through-plane velocity in frame i
        roiDataRaw(:,i) = double(vTemp(roiMask)); %indexed velocities within mask
        meanROI(i) = mean(roiDataRaw(:,i)); %mean velocity in frame i (mm/s)
        flowROI(i) = area.*meanROI(i).*0.001; %flow in frame i (mm^3/s = mL/s)
    end 

    times = double(timeres.*(0:(frames-1))); %original times
%     tq = 0:0.1:size(v,3)-1; %interpolate time dimension
    sampleDensity = 200000;
    tq = linspace(0,frames-1,sampleDensity); %interpolate time dimension
    timesInterp = double(timeres.*tq); %interpolated times

    %Linear interpolation (to get more points on flow curve)
    meanROIfit = interp1(times,meanROI,timesInterp,'linear');
    flowROIfit = interp1(times,flowROI,timesInterp,'linear');

    %Add data to roiStatistics structure
    roiInfo.radius = radius; 
    roiInfo.center = center;
    roiInfo.roiMask = roiMask; 
    roiInfo.roiDataRaw = roiDataRaw;
    roiInfo.Name = ''; %will get changed below
    roiInfo.ROInumber = handles.global.totalROIs;

    %%% Create Interpolated Curve with Gaussian Smoothing           
    meanROIfit = interp1(times,smoothdata(meanROI,'gaussian',3),timesInterp,'cubic');
    flowROIfit = interp1(times,smoothdata(flowROI,'gaussian',3),timesInterp,'cubic');

    roiStatistics.times = timesInterp;
    roiStatistics.meanROI = meanROIfit; 
    roiStatistics.flowROI = flowROIfit; 

    %%% Save all data into handles
    if isstruct(handles.pcDatasets(planeNum).ROIdata)
        handles.pcDatasets(planeNum).ROIinfo(end+1) = roiInfo;
        handles.pcDatasets(planeNum).ROIdata(end+1) = roiStatistics;
    else
        handles.pcDatasets(planeNum).ROIinfo = roiInfo;
        handles.pcDatasets(planeNum).ROIdata = roiStatistics;
    end 

    dataDir = handles.global.homeDir; %directory in which plane data is located
    if dataDir(end)=='\' || dataDir(end)=='/' %kill the slash if it exists
        dataDir(end) = [];
    end 
    if ~exist([dataDir filesep 'ROIimages'],'dir') %if the proposed directory doesn't exist
        mkdir([dataDir filesep 'ROIimages']); %make it
        cd([dataDir filesep 'ROIimages']); %move into it
        frame = getframe(handles.pcPlanePlot); %get a snapshot of the PC plane plot with ROI
        image = frame2im(frame); %make into image
        imwrite(image,[handles.pcDatasets(planeNum).Names '.png']) %write it out as PNG
    else
        cd([dataDir filesep 'ROIimages']); %if ROIimages already exists, move into it
        frame = getframe(handles.pcPlanePlot);
        image = frame2im(frame);
        imwrite(image,[handles.pcDatasets(planeNum).Names '.png'])
    end 
    cd(handles.global.homeDir); %lets go back home
    
    %%% Label each ROI w/ names (helpful because there may be 2 ROIs/plane)
    for i=1:numel(handles.pcDatasets)
        if isstruct(handles.pcDatasets(i).ROIinfo) %if we've made ROI data for this dataset
            if length(handles.pcDatasets(i).ROIinfo)==1
                handles.pcDatasets(i).ROIinfo.Name = handles.pcDatasets(planeNum).Names; %name ROI
            else
                for j=1:length(handles.pcDatasets(i).ROIinfo)
                    name = handles.pcDatasets(i).Names; %get plane name
                    planeName = [name ' ROI ' num2str(j)]; %needed if more than one ROI/plane
                    handles.pcDatasets(i).ROIinfo(j).Name = planeName; %name ROI
                end 
            end 
        end 
    end 

    plotVelocity(handles); %plot flow curves
    
    set(handles.loadCLpush,'Enable','on');
    set(handles.load2DPCbutton,'Enable','on');
    set(handles.pcPlanePopup,'Enable','on'); %make it so we can't select another plane
    set(handles.pcDatasetPopup,'Enable','on'); %make it so we can't select another plane
    set(handles.drawROIbutton,'Enable','on');

    guidata(hObject,handles);
    updatePCImages(handles); %update images (to remove green ROI circle)
    axes(handles.pcPlanePlot); %make sure we're still on PC plot
   

    % --- PLANE SLIDER - CALLBACK
function pcSlider_Callback(hObject, eventdata, handles)
    updatePCImages(handles); %update images if slider is moved

% --- PLANE SLIDER - CREATE FUNCTION
function pcSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
    
% --- MINIMUM CONTRAST VALUE BOX - CALLBACK  
    function minContrastBox_Callback(hObject, eventdata, handles)
    updatePCImages(handles)

% --- MINIMUM CONTRAST VALUE BOX - CREATE FUNCTION   
function minContrastBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- MAXIMUM CONTRAST VALUE BOX - CALLBACK  
function maxContrastBox_Callback(hObject, eventdata, handles)
    updatePCImages(handles)
    
% --- MAXIMUM CONTRAST VALUE BOX - CREATE FUNCTION
function maxContrastBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
    
%%%%%%%%%%%% VELOCITY PLOT %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- VELOCITY PLOT - CREATE FUNCTION
function velocityPlot_CreateFcn(hObject, eventdata, handles)
    

   
    
% --- Executes on button press in completeLoadingROI.
function completeLoadingROI_Callback(hObject, eventdata, handles)
    handles.global.startAnalyzing = 1; %turn on flag to state that we are ready for PWV analysis
    handles.flow = organizeFlowInfo(handles);
    
    % COMPUTE TIME SHIFTS
    if ~isfield(handles.flow,'TTUpstroke')
        flow = computeTTs(handles.flow);
        handles.flow = flow;
    end 
    
    % COMPUTE PWVs
    distance = handles.centerline.PlaneDistances(1);
    TTPoint = [handles.flow(2).TTPoint] - [handles.flow(1).TTPoint]; 
    PWVPoint = distance./TTPoint;
    PWVPoint = PWVPoint(PWVPoint>0);
    PWVPoint = PWVPoint(PWVPoint<40);
    TTFoot = [handles.flow(2).TTFoot] - [handles.flow(1).TTFoot]; 
    PWVFoot = distance./TTFoot;
    PWVFoot = PWVFoot(PWVFoot>0);
    PWVFoot = PWVFoot(PWVFoot<40);
    TTUpstroke = [handles.flow(2).TTUpstroke] - [handles.flow(1).TTUpstroke]; 
    PWVUpstroke = distance./TTUpstroke;
    PWVUpstroke = PWVUpstroke(PWVUpstroke>0);
    PWVUpstroke = PWVUpstroke(PWVUpstroke<40);
    
    PWVP = median(rmoutliers(PWVPoint,'quartiles'));
    PWVF = median(rmoutliers(PWVFoot,'quartiles'));
    PWVU = median(rmoutliers(PWVUpstroke,'quartiles'));
    PWVXC = distance./[handles.flow(2).Xcorr];
    
    numMethods = get(handles.ttpointRadio,'Value') + ...
        get(handles.ttfRadio,'Value') + ...
        get(handles.ttuRadio,'Value') + ...
        get(handles.xcorrRadio,'Value'); %add all PWV buttons turned on
    numCompares = numel(distance); %get number of time shift methods  
    average = zeros(1,numCompares); %initialize average timeshift array
    
    cla(handles.TimeVsDistance,'reset'); %reset PWV plot
    axes(handles.TimeVsDistance); hold on; %force axes to PWV plot
    legendSet = []; %initialize legend cell array
    PWVplots = [];
    iter = 1;
    xlim([0,6]);
    if get(handles.ttpointRadio,'Value')
        PWVplots = [PWVplots PWVPoint];
        legendSet = [legendSet, repmat({'TTPoint'},1,length(PWVPoint))];
        average = average + PWVP; %add all distances for each ROI location
        scatter( iter*ones(size(PWVPoint)),PWVPoint,15,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.03);
        iter = iter+1;
    end 
    if get(handles.ttfRadio,'Value')
        PWVplots = [PWVplots PWVFoot];
        legendSet = [legendSet, repmat({'TTFoot'},1,length(PWVFoot))];
        average = average + PWVF; %keep adding distances
        scatter( iter*ones(size(PWVFoot)),PWVFoot,15,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.03);
        iter = iter+1;
    end 
    if get(handles.ttuRadio,'Value')
        PWVplots = [PWVplots PWVUpstroke];
        legendSet = [legendSet, repmat({'TTUpstroke'},1,length(PWVUpstroke))];
        average = average + PWVU;
        scatter( iter*ones(size(PWVUpstroke)),PWVUpstroke,15,'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.03);
        iter = iter+1;        
    end 
    if get(handles.xcorrRadio,'Value')
        PWVplots = [PWVplots PWVXC];
        legendSet = [legendSet, repmat({'Xcorr'},1,length(PWVXC))];
        average = average + PWVXC;
        iter = iter+1;
    end 
    legendSet{end+1} = 'Avg';
    average = average./numMethods; %get average TT for each ROI
    PWVplots = [PWVplots average];
    boxplot(PWVplots',legendSet')

    if get(handles.ttpointRadio,'Value') %if our ttpoint button is on, plot average ttp
        set(handles.ttpointData,'String',[num2str(round(PWVP,2)) ' m/s']); %write out PWV value in text field
    end 
    if get(handles.ttfRadio,'Value')
        set(handles.ttfData,'String',[num2str(round(PWVF,2)) ' m/s']);
    end 
    if get(handles.ttuRadio,'Value')
        set(handles.ttuData,'String',[num2str(round(PWVU,2)) ' m/s']);
    end 
    if get(handles.xcorrRadio,'Value')
        set(handles.xcorrData,'String',[num2str(round(PWVXC,2)) ' m/s']);
    end 
    
    if numMethods>0 %if we have at least one ttbutton on
        set(handles.averageData,'String',[num2str(round(average,2)) ' m/s']); %set PWV text field
    else %if have not methods selected (all ttbuttons are off)
        set(handles.averageData,'String','0 m/s'); %set PWV text field to 0 m/s
    end
    
    set(handles.ttpointRadio,'Enable','on');
    set(handles.ttuRadio,'Enable','on');
    set(handles.ttfRadio,'Enable','on');
    set(handles.xcorrRadio,'Enable','on');
    set(handles.exportAnalysisButton,'Enable','on');
    guidata(hObject, handles);
    
    

%%%%%%%%%%%% PWV PLOT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- PWV PLOT - CREATE FUNCTION
function TimeVsDistance_CreateFcn(hObject, eventdata, handles)


% --- EXPORT ANALYSIS - CALLBACK
function exportAnalysisButton_Callback(hObject, eventdata, handles)
    flow = computeTTs(handles.flow); %recalculate flow struct 
    numROIs = numel(flow);

    Distance = handles.centerline.PlaneDistances(1);
    TTPoint = [handles.flow(2).TTPoint] - [handles.flow(1).TTPoint]; 
    PWVPoint = Distance./TTPoint;
    PWVPoint = PWVPoint(PWVPoint>0);
    PWVPoint = PWVPoint(PWVPoint<40);
    TTFoot = [handles.flow(2).TTFoot] - [handles.flow(1).TTFoot]; 
    PWVFoot = Distance./TTFoot;
    PWVFoot = PWVFoot(PWVFoot>0);
    PWVFoot = PWVFoot(PWVFoot<40);
    TTUpstroke = [handles.flow(2).TTUpstroke] - [handles.flow(1).TTUpstroke]; 
    PWVUpstroke = Distance./TTUpstroke;
    PWVUpstroke = PWVUpstroke(PWVUpstroke>0);
    PWVUpstroke = PWVUpstroke(PWVUpstroke<40);
    Xcorr = [handles.flow(2).Xcorr];
    
    
    TTP = median(rmoutliers(TTPoint,'quartiles'));
    TTF = median(rmoutliers(TTFoot,'quartiles'));
    TTU = median(rmoutliers(TTUpstroke,'quartiles'));
    TTAverage = mean([TTP; TTF; TTU; Xcorr],1); %get average time shift
    
    PWVP = median(rmoutliers(PWVPoint,'quartiles'));
    PWVF = median(rmoutliers(PWVFoot,'quartiles'));
    PWVU = median(rmoutliers(PWVUpstroke,'quartiles'));
    PWVXC = Distance./Xcorr;
    PWVAv = Distance./TTAverage;

    % Make first column for excel file
    PLANES = [flow(1).Name ' --> ' flow(2).Name]; %get names for Excel

    %PLANES = PLANES'; %Needed for excel, won't save name if ' in table call
    Distance = Distance';
    TTPoint = TTPoint';
    TTFoot = TTFoot';
    TTUpstroke = TTUpstroke';
    Xcorr = Xcorr';
    TTAverage = TTAverage';

    % Make table for writing excel file
    pwvTable = table({PLANES},Distance,TTP,TTF,TTU,Xcorr,TTAverage,PWVP,PWVF,PWVU,PWVXC,PWVAv);
    baseDir = handles.global.homeDir; %rejoin string to get name of folder one up from plane data
    date = datestr(now); %get current date/time
    chopDate = [date(1:2) '-' date(4:6) '-' date(10:11) '-' date(13:14) date(16:17)]; %chop date up
    if ~exist([baseDir filesep 'DataAnalysis'],'dir') %if directory doesn't exist
        mkdir([baseDir filesep 'DataAnalysis']); %make it
        cd([baseDir filesep 'DataAnalysis']); %go to it
        writetable(pwvTable,['Summary_' chopDate '.xlsx'],'FileType','spreadsheet'); %write excel sheet
        flow = handles.flow; %make variables for saving
        save('flow.mat','flow')
        save('pwvTable.mat','pwvTable')
    else %or if the directory already exists
        cd([baseDir filesep 'DataAnalysis']); %go to it
        writetable(pwvTable,['Summary_' chopDate '.xlsx'],'FileType','spreadsheet');
        flow = handles.flow;
        save('flow.mat','flow')
        save('pwvTable.mat','pwvTable')
    end 
    cd(handles.global.homeDir); %go back home  
    clear PLANES %need to do this because PLANES will keep getting transposed
    
    cd([baseDir filesep 'DataAnalysis']); %go to it
    frame = getframe(handles.TimeVsDistance); %get snapshot of PWV plot 
    imwrite(frame2im(frame),'PWVanalysisPlot.png'); %write out to PNG
    pcDatasets = handles.pcDatasets;
    save('pcDatasets.mat','pcDatasets');
    cd(handles.global.homeDir); %go back home  
    mkdir('PWV_RT_Analysis');
    movefile('DataAnalysis','PWV_RT_Analysis');
    movefile('ROIimages','PWV_RT_Analysis');
    set(handles.exportDone,'String','Export Completed!');

    guidata(hObject, handles);
    
   

% --- TTPoint READOUT - CREATE FUNCTION
function ttpointData_CreateFcn(hObject, eventdata, handles)
% --- TTPoint RADIO - CALLBACK
function ttpointRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttpointRadio,'Value') %if we're turned off
        set(handles.ttpointData,'String',' '); %don't display PWV
    end 
    if handles.global.startAnalyzing %if we're analyzing PWVs
        completeLoadingROI_Callback(hObject, eventdata, handles); %reanalyze without TTpoint 
    end 

% --- TTUpstroke READOUT - CREATE FUNCTION
function ttuData_CreateFcn(hObject, eventdata, handles)
% --- TTUpstroke RADIO - CALLBACK
function ttuRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttuRadio,'Value')
        set(handles.ttuData,'String',' ');
    end 
    if handles.global.startAnalyzing
        completeLoadingROI_Callback(hObject, eventdata, handles);
    end 

% --- TTFoot READOUT - CREATE FUNCTION
function ttfData_CreateFcn(hObject, eventdata, handles)
% --- TTFoot RADIO - CALLBACK
function ttfRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttfRadio,'Value')
        set(handles.ttfData,'String',' ');
    end 
    if handles.global.startAnalyzing
        completeLoadingROI_Callback(hObject, eventdata, handles);
    end 

% --- Xcorr READOUT - CREATE FUNCTION
function xcorrData_CreateFcn(hObject, eventdata, handles)
% --- Xcorr RADIO - CALLBACK
function xcorrRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.xcorrRadio,'Value')
        set(handles.xcorrData,'String',' ');
    end 
    if handles.global.startAnalyzing
        completeLoadingROI_Callback(hObject, eventdata, handles);
    end 

% --- AVERAGE READOUT - CREATE FUNCTION
function averageData_CreateFcn(hObject, eventdata, handles)





%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% --- Update Images in PLANE PLOT
function updatePCImages(handles)
    axes(handles.pcPlanePlot); %force axes to PC plot
    
    planeNum = get(handles.pcPlanePopup,'Value'); % get current plane index (eg AAo)
    datasetNum = get(handles.pcDatasetPopup,'Value'); %get current dataset (eg MAG)
    
    dataset = handles.pcDatasets;
    if isstruct(dataset(planeNum).Images) %if Data is a structure
        imageSet = struct2cell(dataset(planeNum).Images); %convert from struct to cell
    else  %or if we already have an array
        imageSet = dataset(planeNum).Images; %just change the name
    end 

    if iscell(imageSet) %if our set of images are contained in a cell
        images = imageSet(datasetNum); %pull images for current dataset
        images = cell2mat(images); %turn to matrix
    else 
        images = imageSet; %otherwise, do nothing
    end 
    
    if ndims(images)<3 %if we are dealing with time-averaged images (ndim=2)
        maxSize = max(size(images,1),size(images,2));
        steps = [1 maxSize]; %set our slider as wide as possible so we can't slide
        set(handles.pcSlider,'SliderStep', steps);
        slice = images; %change name for consistency (see below)
    else %if we are dealing with time-resolved images (ndim=3)
        dim3size = size(images,3); %get the size of the third dimension
        steps = [1/(dim3size-1) 10/(dim3size-1)]; %set so one 'slide' moves to the next slice exactly
        set(handles.pcSlider,'SliderStep', steps);
        sliceNum = 1+round( get(handles.pcSlider,'Value').*(dim3size-1) ); %get slice number from slider
        slice = images(:,:,sliceNum); %pull slice from images
    end 
    
    minc = str2double(get(handles.minContrastBox,'String'));
    maxc = str2double(get(handles.maxContrastBox,'String'));

    if ~isempty(handles.pcDatasets(planeNum).ROI) %if we have an ROI placed
        hold on %make it so we can slide while keeping the ROI on the figure
        imshow(rescale(slice),[minc maxc]);
    else 
        cla(handles.pcPlanePlot,'reset') %otherwise, reset the plot
        imshow(rescale(slice),[minc maxc]) %then show the image
    end 
    
    
    
% --- "Time to" calculations (TTPoint, TTUpstroke, TTFoot, Xcorr)
function flow = computeTTs(flow)
    v1 = [flow(1).Gaussian.flowROI];
    v2 = -1*[flow(2).Gaussian.flowROI];
    times = [flow(1).Gaussian.times];
    timeres = times(2)-times(1);

    minRRdist = 350;
    %60000*tempRes/maxHR = minRRdist;
    M1 = max(v1);
    [~,p1] = findpeaks(v1,'MinPeakHeight',M1/2,'MinPeakDistance',minRRdist);
    peaks1 = times(p1);
    TTPeak(1,:) = peaks1;

    for k=1:length(p1)
        diff = 1;
        it = 0;
        while diff>0
            diff = v1(p1(k)-it) - v1(p1(k)-it-1);
            it = it+1;
        end 
        min1(k) = p1(k)-it+1;
        t0 = times(min1(k):p1(k));
        upslope = v1(min1(k):p1(k));

        % TTPoint - time to point calculation
        [~,idx] = min(abs(rescale(upslope)-0.5));
        TTPoint(1,k) = t0(idx);

        % TTUpstroke - time to upstroke calculation
        [~,~,t1] = sigFitRT(upslope,t0);
        TTUpstroke(1,k) = t1;

        % TTFoot - time to foot calculation
        [~,leftIdx] = min(abs(rescale(upslope)-0.2));
        [~,rightIdx] = min(abs(rescale(upslope)-0.8));
        vSeg = upslope(leftIdx:rightIdx);
        tSeg = t0(leftIdx:rightIdx);
        m = polyfit(tSeg,vSeg,1);
        TTFoot(1,k) = -(m(2)/m(1));
    end 
    flow(1).TTPeak = TTPeak(1,:);
    flow(1).TTPoint = TTPoint(1,:);
    flow(1).TTUpstroke = TTUpstroke(1,:);
    flow(1).TTFoot = TTFoot(1,:);

    M2 = max(v2);
    [~,p2] = findpeaks(v2,'MinPeakHeight',M2/2,'MinPeakDistance',minRRdist);
    peaks2 = times(p2);
    TTPeak(2,:) = peaks2;

    for k=1:length(p2)
        diff = 1;
        it = 0;
        while diff>0
            diff = v2(p2(k)-it) - v2(p2(k)-it-1);
            it = it+1;
        end 
        min2(k) = p2(k)-it+1;
        t0 = times(min2(k):p2(k));
        upslope = v2(min2(k):p2(k));

        % TTPoint - time to point calculation
        [~,idx] = min(abs(rescale(upslope)-0.5));
        TTPoint(2,k) = t0(idx);

        % TTUpstroke - time to upstroke calculation
        [~,~,t1] = sigFitRT(upslope,t0);
        TTUpstroke(2,k) = t1;

        % TTFoot - time to foot calculation
        [~,leftIdx] = min(abs(rescale(upslope)-0.2));
        [~,rightIdx] = min(abs(rescale(upslope)-0.8));
        vSeg = upslope(leftIdx:rightIdx);
        tSeg = t0(leftIdx:rightIdx);
        m = polyfit(tSeg,vSeg,1);
        TTFoot(2,k) = -(m(2)/m(1));
    end 
    flow(2).TTPeak = TTPeak(2,:);
    flow(2).TTPoint = TTPoint(2,:);
    flow(2).TTUpstroke = TTUpstroke(2,:);
    flow(2).TTFoot = TTFoot(2,:);
    
    % XCorr - cross correlation calculation
    [c,lags] = xcorr(v2,v1);
    [~,maxC] = max(c); %get index of max Xcorr value
    shift = lags(maxC); %find time lag of Xcorr peak
    flow(1).Xcorr = 0;
    flow(2).Xcorr = shift*timeres;


 
% --- Turn PolyLine into SplineLine    
function Y = interppolygon(X,N)
    if nargin < 2 || N < 2
        N = 2; %if only one arg or too small N, just assume 2
    end
    nDim = size(X,2); %should be 2
    dx = 0;
    
    for dim = 1:nDim
        dx = dx + diff(X(:,dim)).^2 ; %get sum of squares in each dim
    end
    
    lengthBetweenPoints = sqrt(dx); %now get distance
    lengthLine = sum(lengthBetweenPoints);
    origMetric = [0; cumsum(lengthBetweenPoints/lengthLine)];
    
    interpMetric = (0:(1/(N-1)):1)';
    Y = interp1(origMetric,X,interpMetric,'makima'); %makima seems to work well
    %Y = csaps([0 times times(end)+times(1)],[0 meanROI 0],0.0001,timesInterp);
    

    
% --- Plot Velocities    
function plotVelocity(handles)
    cla(handles.velocityPlot,'reset'); %reset axes
    axes(handles.velocityPlot); %make sure we plot on the right axis
    
    times = handles.pcDatasets(1).ROIdata(1).times;
    plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5); %line of y=0 (for visual reference)
    % Note that above we assume same time scale for each plane (MR scan)
    xlim([min(times) max(times)]);
    legendSet = {'Baseline'}; %add baseline to legend names
    for i=1:numel(handles.pcDatasets) %for all planes
        if isstruct(handles.pcDatasets(i).ROIdata) %if we've made ROI data for this dataset
            for j=1:length(handles.pcDatasets(i).ROIdata) %for each ROI
                legendSet{end+1} = handles.pcDatasets(i).ROIinfo(j).Name; %add name of ROI to list  
                times = handles.pcDatasets(i).ROIdata(j).times;
                velocity = handles.pcDatasets(i).ROIdata(j).meanROI; %grab mean velocity
                if mean(velocity)<0 %if we are mainly negative velocities (as in descending aorta)
                    hold on; plot(times,-1*velocity); %else, plot inverted velocity
                else %otherwise, don't invert velocity (as in ascending aorta)
                    hold on; plot(times,velocity);
                end 
            end 
        end 
    end 
    legend(legendSet); hold off
    xlabel('Time (ms)'); ylabel('Mean Velocity in ROI (mm/s)'); %set axes labels
    %xlim([times(1),times(end)]); %chop limits to make curve full width
    

    
% --- Sigmoid Fit Function
function [sigmoid,t,t1] = sigFitRT(upslope,t0)
%%% See the following article by Anas Dogui in JMRI:
% Measurement of Aortic Arch Pulse Wave Velocity in Cardiovascular MR:
% Comparison of Transit Time Estimators and Description of a New Approach

    t = linspace(t0(1),t0(end),1000); %interpolate even more
    upslope = interp1(t0,upslope,t); %interpolate upslope
    upslope = rescale(upslope); %normalize from 0 to 1
    dt = t(2)-t(1); %new temporal resolution (=0.1)
    midpoint = round(length(upslope)/2);
    
    % c1 = b, c2 = a, c3 = x0, c4 = dx
    % Note that we could assume the equation e^t/(1+e^(t-t0)) since c1=1 and
    % c2=0. However, will keep the same as the Dogui paper.
    sigmoidModel = @(c) c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ) - upslope;
    c0 = [0,1,t(midpoint),dt/2]; %initial params for upslope region
    opts = optimset('Display', 'off'); %turn off display output
    c = lsqnonlin(sigmoidModel,c0,[],[],opts); %get nonlinear LSQ solution
    sigmoid = c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ); %calculate our sigmoid fit with our new params

    dy = diff(sigmoid,1);
    dy(end) = [];
    dx = diff(t,1);
    dx(end) = [];
    ddy = diff(sigmoid,2);
    curvature = ddy.*dx./(dx.^2 + dy.^2).^(3/2);
    [~,tIdx] = max(curvature);
    t1 = t(tIdx);

    
    
% --- Condense and organize flow data obtained from ROIs
function flow = organizeFlowInfo(handles)
% This function is designed to pull apart handles.pcDatasets. It is
% difficult to do analysis on the handles structure because some slices
% have 2 ROIs. It is much easier to have a structure that pulls out each
% ROI. It only adds a bit of memory since we aren't saving the raw images.
    count = 1; %overall iterator
    for i=1:numel(handles.pcDatasets) %for each PC dataset
        if isstruct(handles.pcDatasets(i).ROIdata) %do we have data?
            if length(handles.pcDatasets(i).ROIdata)==1 %if we just have one ROI
                flow(count).Name = handles.pcDatasets(i).ROIinfo.Name; %pull only relevant info for PWV calcs
                flow(count).ROIinfo = handles.pcDatasets(i).ROIinfo;
                flow(count).HeaderInfo = handles.pcDatasets(i).Info;
                flow(count).Gaussian = handles.pcDatasets(i).ROIdata;
                flow(count).pcDatasetREF = [i,1];
                count = count+1;
            else 
                for j=1:length(handles.pcDatasets(i).ROIdata) %if we have more than one ROI
                    flow(count).Name = handles.pcDatasets(i).ROIinfo(j).Name; %parse into individual ROIs in flow struct
                    flow(count).ROIinfo = handles.pcDatasets(i).ROIinfo(j);
                    flow(count).HeaderInfo = handles.pcDatasets(i).Info;
                    flow(count).Gaussian = handles.pcDatasets(i).ROIdata(j);
                    flow(count).pcDatasetREF = [i,j];
                    count = count+1;
                end 
            end 
        end 
    end 

 function saveTTplots(handles,flow)   
    times = handles.pcDatasets(1).ROIdata(1).times;
    figure('units','normalized','outerposition',[0 0 1 1]); 
    plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5); %line of y=0 (for visual reference)
    % Note that above we assume same time scale for each plane (MR scan)
    xlim([min(times) max(times)]);
    legendSet = {'Baseline'}; %add baseline to legend names
    count = 1;
    for i=1:numel(handles.pcDatasets) %for all planes
        if isstruct(handles.pcDatasets(i).ROIdata) %if we've made ROI data for this dataset
            for j=1:length(handles.pcDatasets(i).ROIdata) %for each ROI
                legendSet{end+1} = handles.pcDatasets(i).ROIinfo(j).Name; %add name of ROI to list
                vTemp = handles.pcDatasets(i).ROIdata(j).meanROI;
                
                if mean(vTemp)<0
                    vTemp = -1*vTemp;
                end
                
                hold on; plot(times,vTemp);
                velocity(count,:) = vTemp;
                count = count+1;
            end 
        end 
    end 
    legend(legendSet,'Location','northeastoutside','AutoUpdate','off');
    xlabel('Time (ms)'); 
    mx = max(velocity(:));
    mn = min(velocity(:));
    ylim([(mn - mn*0.05) (mx + mx*0.05)]);
    ylabel('Mean Velocity in ROI (mm/s)'); %set axes labels
    ttp = [flow.TTPoint];
    for i=1:length(ttp)
        [~,idx] = min( abs(times-ttp(i)) );
        h(i) = scatter(times(idx),velocity(i,idx),'k','filled');
        h2(i) = scatter(times(idx),0,'k');
        h3(i) = plot([times(idx) times(idx)],[velocity(i,idx),0],':k');
    end 
    title('Time-To-Point');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'TTPoint.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off'); set(h3,'Visible','off')
    
    ttf = [flow.TTFoot];
    for i=1:length(ttf)
        [~,idx] = min( abs(times-ttf(i)) );
        [scale,idx2] = max(velocity(i,:));
        P1 = flow(i).P1;
        x = times(idx:idx2);
        y = scale*(P1(1)*x + P1(2));
        h(i) = plot(x,y,'k');
        h2(i) = scatter(times(idx),0,'k');
    end 
    title('Time-To-Foot');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'TTFoot.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off');
        
    ttu = [flow.TTUpstroke];
    for i=1:length(ttu)
        [~,idx] = min( abs(times-ttu(i)) );
        [scale,~] = max(velocity(i,:));
        h(i) = scatter(times(idx),velocity(i,idx),'k','filled');
        h2(i) = scatter(times(idx),0,'k');
        h3(i) = plot([times(idx) times(idx)],[velocity(i,idx),0],':k');
        h4(i) = plot(flow(i).SigmoidTimes,scale*flow(i).SigmoidFit,'.k');
    end 
    title('Time-To-Upstroke');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'TTUpstroke.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off'); set(h3,'Visible','off'); set(h4,'Visible','off');
    
    xcorr = [flow.Xcorr];
    for i=1:length(xcorr)
        [~,idx] = min( abs(times-xcorr(i)) );
        h(i) = scatter(times(idx),velocity(i,idx),'k','filled');
        h2(i) = scatter(times(idx),0,'k');
        h3(i) = plot([times(idx) times(idx)],[velocity(i,idx),0],':k');
    end 
    title('Cross Correlation Time Lag');
    frame = getframe(gcf);
    imwrite(frame2im(frame),[pwd filesep 'Xcorr.png']);
    
    set(h,'Visible','off'); set(h2,'Visible','off'); set(h3,'Visible','off')
    
    hold off; close(gcf);   
    
    
    
% Load Dat files
function v = load_dat(name, res)
    [fid,errmsg]= fopen(name,'r');
    if fid < 0  %if name does not exist in directory
        set(handles.MessageBar,'String',['Error Opening Data : ',errmsg]);
    end

    % Reads in as short, reshapes by image res.
    v = reshape(fread(fid,'short=>single'),res);
    fclose(fid);
