%---------------------- Basic Figure ---------------------%
fh = figure('Units', 'Pixels', 'OuterPosition', [100 100 1500 900], 'Toolbar', 'none', 'Menu', 'none','Color','w','Tag','MainWindow');


%---------------- Axes and Image Controls ---------------%
H.axes1 = axes('Units','pixels','Tag','axes1','Position',[600 50 850 750],'Box','on','XTick',[],'YTick',[]); pause(H.delay);        
H.image = imagesc(zeros(100,100),[0 1]);
colormap(CubeHelix(256,0.5,-1.5,1.2,1.0));
set(gca,'Box','on','XTick',[],'YTick',[]);

    %----- Sliders -----%
    H.slider.frame = uicontrol( 'Style', 'slider','Position',[600 30 850 19],'Background', 'w', 'Callback', @FrameSlider); pause(H.delay);            
    H.slider.min = uicontrol( 'Style', 'slider', 'Parent', fh, 'Background', 'w', 'Position',[1451 50 19 750], 'Callback', @ColorRange); pause(H.delay);       
    H.slider.max = uicontrol( 'Style', 'slider', 'Parent', fh, 'Position',[1471 50 19 750],'Background', 'w','Callback', @ColorRange); pause(.5);       

    %----- Radio Buttons -----%
    H.rdio.stack = uicontrol( 'Style', 'radiobutton', 'Parent', fh,'Callback', @ShowWhat); pause(H.delay);
    set(H.rdio.stack,'Position',[600 801 100 25], 'Background', 'w', 'String','Show Stack');
    H.rdio.zproj = uicontrol( 'Style', 'radiobutton', 'Parent', fh,'Callback', @ShowWhat); pause(H.delay);
    set(H.rdio.zproj,'Position',[700 801 125 25], 'Background', 'w', 'String','Show Projection');
    H.rdio.roi = uicontrol( 'Style', 'radiobutton', 'Parent', fh,'Callback', @ShowWhat); pause(H.delay);
    set(H.rdio.roi,'Position',[825 801 100 25], 'Background', 'w', 'String','Show ROIs');
    
    %----- Popup Menus -----%
    H.popup.roitype = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'Position',[1350 801 100 30], 'Background', 'w', 'String',{'Scatter';'T-Scatter';'Patch';'Outline'}); pause(H.delay);
    H.popup.zprojtype = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'Position',[1225 801 100 30], 'Background', 'w', 'String',{'Mean';'STD';'Max';'CC';'MAD';'Kurtosis';'Mean AND STD';'Mean OR STD'}); pause(H.delay);

    
%--------------- Load and Save Controls ---------------%

    %----- Buttons -----%
    H.button.load = uicontrol( 'Style', 'pushbutton', 'Parent', fh,'Callback', @LoadStack); pause(H.delay);
    set(H.button.load,'Position',[25 835 100 30], 'Background', 'w', 'String','Load');
    H.button.save = uicontrol( 'Style', 'pushbutton', 'Parent', fh,'Callback', @SaveStack); pause(H.delay);
    set(H.button.save,'Position',[130 835 100 30], 'Background', 'w', 'String','Save');

    %----- Path Edit String -----%
    H.edit.filepath = uicontrol( 'Style', 'edit', 'Parent', fh); pause(H.delay);
    set(H.edit.filepath,'Position',[120 805 450 25], 'Background', 'w');            

    %----- Drop Menu -----%
    H.pop.stacktype = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'String',{'Raw';'Registered';'Processed'}); pause(H.delay);
    set(H.pop.stacktype,'Position',[25 805 90 25], 'Background', 'w'); 

    
%--------------- Stack Registration ---------------%
H.button.reg = uicontrol( 'Style', 'pushbutton', 'Parent', fh,'Callback', @Registration,'Position',[25 750 100 30], 'Background', 'w', 'String','Register'); pause(H.delay);


%--------------- Filter ---------------%
H.button.filt = uicontrol( 'Style', 'pushbutton', 'Parent', fh,'Callback', @Filter,'Position',[25 700 100 30], 'Background', 'w', 'String','Filter'); pause(H.delay);
H.popup.filthow = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'Position',[130 700 100 30], 'Background', 'w', 'String',{'Kalman';'Median';'Wiener';'Original'}); pause(H.delay);
H.popup.filtwhat = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'Position',[235 700 100 30], 'Background', 'w', 'String',{'Stack';'Proj'}); pause(H.delay);
H.popup.wproj = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'Position',[340 700 100 30], 'Background', 'w', 'String',{'Mean';'STD';'CC'}); pause(H.delay);


%--------------- Get ROIs and Segmentation ---------------%
    
    %---------- Automated Routines ----------%
    H.text.autoroilabel = uicontrol( 'Style', 'text', 'Parent', fh,'Position',[25 650 100 25], 'Background', 'w', 'String','Auto ROIs','FontWeight','bold'); pause(H.delay);
    H.button.autoroi = uicontrol( 'Style', 'pushbutton', 'Parent', fh,'Callback', @AutoRoi,'Position',[25 625 100 30], 'Background', 'w', 'String','Get ROIs'); pause(H.delay);
    H.popup.roihow = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'Position',[130 625 100 30], 'Background', 'w', 'String',{'Laplace';'Threshold';'Correlation';'ICA';'NCuts'}); pause(H.delay);
    H.popup.roiwhat = uicontrol( 'Style', 'popupmenu', 'Parent', fh,'Position',[235 625 100 30], 'Background', 'w', 'String',{'Mean';'STD';'CC'}); pause(H.delay);
    H.edit.roiparam = uicontrol( 'Style', 'edit', 'Parent', fh,'Position',[340 625 250 30], 'Background', 'w','String',''); pause(H.delay);

    %---------- Manual Routines ----------%
    H.text.autoroilabel = uicontrol( 'Style', 'text', 'Parent', fh,'Position',[25 575 100 25], 'Background', 'w', 'String','Manual ROIs','FontWeight','bold'); pause(H.delay);
    H.button.initroi = uicontrol( 'Style', 'pushbutton', 'Parent', fh,'Callback', @InitRoi,'Position',[25 550 100 30], 'Background', 'w', 'String','Get ROIs'); pause(H.delay);
    H.rdio.addroi = uicontrol( 'Style', 'radiobutton', 'Parent', fh,'CallBack',@AddRoi,'Position',[25 522 100 30], 'Background', 'w', 'String','Add ROI'); pause(H.delay);
    H.rdio.removeroi = uicontrol( 'Style', 'radiobutton', 'Parent', fh,'Position',[25 500 100 30], 'Background', 'w', 'String','Remove ROI'); pause(H.delay);
    

        
%--------------- Notifications ---------------%        
H.edit.activity = uicontrol( 'Style', 'edit', 'Parent', fh,'Position',[600 835 600 25], 'Background', 'w','String',''); pause(H.delay);

        
        
        