function ImageViewerK()
global H
H = [];
%-------------------- Draw GUI & Components --------------------%
H.delay = .02;
MakeGUI;    
                           
%-------------------- Initiate Variables --------------------%
InitiateVariables;

            
%|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||%
%-------------------------- Sub Functions --------------------------------%
%|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||%

%--------------- Load and Save Functions ---------------%
    
    %---------- Load Stack ----------%
    function LoadStack(~,~)
    global H    
        %---------- Get FilePath ----------%
        H.FilePath = uigetdir(H.filepath,'Get Directory for Data');

        %---------- Load Image Stack ----------%
        [H.ZStack.c, H.NFrames]=readtiff2(H.FilePath);      

        %---------- Make Projections ----------%
        H.zproj.m = mean(double(H.ZStack.c),3);
        H.zproj.s = std(double(H.ZStack.c),1,3);
        H.zproj.k = kurtosis(double(H.ZStack.c),0,3);
        H.zproj.max = max(double(H.ZStack.c),[],3);
        H.zproj.mad = mad(double(H.ZStack.c),0,3);             

        %---------- Set Some Variables ----------%
        H.Max = max(H.ZStack.c(:));
        H.Min = min(H.ZStack.c(:));   
        H.ZStack.org = H.ZStack.c;   

        %---------- Display First Image of Stack ----------%
        axes(H.axes1)
        H.image = imagesc(H.ZStack.c(:,:,1),[H.Min H.Max]); colormap(CubeHelix(256,0.5,-1.5,1.2,1.0))
        set(gca,'Box','on','XTick',[],'YTick',[])

        %---------- Send Info to GUI ----------%
        set(H.slider.frame,'Max',H.NFrames,'Min',1,'SliderStep',[1/H.NFrames 10/H.NFrames],'Value',1);
        set(H.slider.min,'Min',0,'Max',H.Max,'Value',H.Min);
        set(H.slider.max,'Min',0,'Max',H.Max,'Value',H.Max);
        set(H.rdio.zproj,'Value',0);
        set(H.rdio.roi,'Value',0);
        set(H.rdio.stack,'Value',1);
        H.ShowWhat = 'Stack';
                                   

%--------------- Image Display and Manipulation ---------------%

    %---------- Show Which Image ----------%
    function ShowWhat(hObj,~)
    global H
    axes(H.axes1);

    %---------- Get Event ----------%
    W = get(hObj,'String');
    switch W
        case 'Show Stack'
            set(H.rdio.stack,'Value',1);
            set(H.rdio.zproj,'Value',0);
            set(H.rdio.roi,'Value',0);
            H.ShowWhat = 'Stack';
            H.Max = max(H.ZStack.c(:));
            H.Min = min(H.ZStack.c(:));            
            H.image = imagesc(H.ZStack.c(:,:,1),[H.Min H.Max]); 
            set(H.slider.frame,'Max',H.NFrames,'Min',1,'SliderStep',[1/H.NFrames 10/H.NFrames],'Value',1);
        case 'Show Projection'
            set(H.rdio.stack,'Value',0);
            set(H.rdio.zproj,'Value',1);
            set(H.rdio.roi,'Value',0);
            H.ShowWhat = 'Proj';
            wproj = get(H.popup.zprojtype,'String'); wproj = wproj{get(H.popup.zprojtype,'Value')};
            switch wproj
                case 'Mean'
                    Proj = H.zproj.m;
                case 'STD'
                    Proj = H.zproj.s;
                case 'Kurtosis'
                    Proj = H.zproj.s;
                case 'MAD'
                    Proj = H.zproj.mad;
                case 'Max'
                    Proj = H.zproj.max;
                case 'Mean AND STD'
                    Proj = H.zproj.m .* H.zproj.s;   
                case 'CC'
                    if ~isfield(H.zproj,'c')
                        pp = gcp('nocreate');
                        if isempty(pp)
                            parpool(H.npools); H.poolobj = gcp('nocreate');
                        end
                        H.zproj.c = CrossCorrImageKP(double(H.ZStack.c),H.npools,'CC','yes','no');
                    end
                    Proj = H.zproj.c;
                case  'Mean AND STD'
                    disp('Program it')
            end
            Proj = (Proj - min(Proj(:)))/(max(Proj(:)) - min(Proj(:)));
            H.Max = max(Proj(:));
            H.Min = min(Proj(:));
            H.image = imagesc(Proj,[H.Min H.Max]);
            set(H.slider.frame,'Max',2,'Min',1,'SliderStep',[1 1],'Value',1);
        case 'Show ROIs'
            wproj = get(H.popup.roiwhat,'String'); wproj = wproj{get(H.popup.roiwhat,'Value')};
            switch wproj
                case 'Mean'
                    Proj = H.zproj.m;
                case 'STD'
                    Proj = H.zproj.s;
                case 'CC'
                    if ~isfield(H.zproj,'c')
                        pp = gcp('nocreate');
                        if isempty(pp)
                            parpool(H.npools); H.poolobj = gcp('nocreate');
                        end
                        H.zproj.c = CrossCorrImageKP(double(H.ZStack.c),H.npools,'CC','yes','no');
                    end
                    Proj = H.zproj.c;
            end
            Proj = (Proj - min(Proj(:)))/(max(Proj(:)) - min(Proj(:)));
            set(H.rdio.stack,'Value',0);
            set(H.rdio.zproj,'Value',0);
            set(H.rdio.roi,'Value',1);
            H.Max = max(Proj(:));
            H.Min = min(Proj(:));
            roitype = get(H.popup.roitype,'String'); roitype = roitype{get(H.popup.roitype,'Value')};
            axes(H.axes1)
            H.image = imagesc(Proj); hold on
            colormap(CubeHelix(256,0.5,-1.5,1.2,1.0))
            set(gca,'Box','on','XTick',[],'YTick',[]);
            for i = 1:H.roi.con.NumObjects            
                switch roitype
                    case 'Scatter'                    
                        scatter(H.roi.reg(i).PixelList(:,1),H.roi.reg(i).PixelList(:,2),50,H.roiCM(i,:),'fill','MarkerEdgeColor','none');
                    case 'T-Scatter'
                        scatter_patches(H.roi.reg(i).PixelList(:,1),H.roi.reg(i).PixelList(:,2),50,H.roiCM(i,:),'FaceAlpha',0.1,'EdgeColor','none');
                    case 'Patch'                    
                        patch(H.roi.reg(i).ConvexHull(:,1),H.roi.reg(i).ConvexHull(:,2),'-','FaceColor',H.roiCM(i,:),'FaceAlpha',0.5,'EdgeColor','none');                    
                    case 'Outline'
    %                     contour(H.roi.reg(i).Image,1,'Color',H.roiCM(i,:));                    
                        plot(H.roi.reg(i).ConvexHull(:,1),H.roi.reg(i).ConvexHull(:,2),'-','Color',H.roiCM(i,:),'LineWidth',2);                    
                end            
            end 
            hold off
            set(H.slider.frame,'Max',2,'Min',1,'SliderStep',[1 1],'Value',1);
        otherwise
    end
    colormap(CubeHelix(256,0.5,-1.5,1.2,1.0));
    set(gca,'Box','on','XTick',[],'YTick',[]);

    %---------- Clean up GUI ----------%
    set(H.slider.min,'Min',0,'Max',H.Max,'Value',H.Min);
    set(H.slider.max,'Min',0,'Max',H.Max,'Value',H.Max);    
    
    %---------- Frame Slider ----------%
    function FrameSlider(~,~)
    global H
    
        frame = round(get(H.slider.frame,'Value'));
        if frame < 1 
            frame = 1;
        elseif frame > H.NFrames
            frame = H.NFrames;
        end                
        set(H.image,'CData',H.ZStack.c(:,:,frame));
        colormap(CubeHelix(256,0.5,-1.5,1.2,1.0));
        set(gca,'Box','on','XTick',[],'YTick',[]);
        
    %---------- Set Color Limits ----------%
    function ColorRange(~,~)
    global H
    
        M(1) = get(H.slider.min,'Value');
        M(2) = get(H.slider.max,'Value');        
        switch H.ShowWhat            
            case 'Stack'
                frame = round(get(H.slider.frame,'Value'));
                H.image = imagesc(H.ZStack.c(:,:,frame),M);  
            case 'Proj'
                H.image = imagesc(H.zproj.m,M);
            case 'ROI'
                H.image = imagesc(H.ROI,M);
        end
        

%--------------- Stack Filter and Registration ---------------%

    %---------- Registration ----------%
    function Registration(~,~)
    global H
        
        %---------- Register Stack ----------%
        p = [1 10 10];        
        target = H.zproj.m;
        defaultopts = {'nStripes',p(1),'padLen',p(2)};
        options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
        nstacks = H.npools;
        stack = double(H.ZStack.c);
        [~, ~, nz] = size(stack);
        nframes = zeros(1,nstacks);
        if rem(nz,nstacks) == 0
            nframes(:) = nz/nstacks;
        elseif rem(nz,nstacks) == 1 
            nframes(end) = ceil(nz/nstacks);
            nframes(1:end-1) = (nz - nframes(end))/(nstacks - 1); 
        else
            nframes(1:end-1) = floor(nz/nstacks);
            nframes(end) = nz - sum(nframes(1:end-1));
        end
        smallStacks = cell(nstacks,1);
        index1 = 1;
        index2 = nframes(1);
        smallStacks{1} = stack(:,:,index1:index2);
        for i = 2:nstacks                
            index1 = (i-1)*nframes(i-1)+1;
            index2 = sum(nframes(1:i));
            smallStacks{i,1} = stack(:,:,index1:index2);
        end
        parfor i = 1:nstacks
            smallStacks{i,1} = stackTurboReg2(smallStacks{i,1},target,options.nStripes,options.padLen);
        end
        bigStack = [];
       for i = 1:nstacks
            bigStack = cat(3,bigStack,smallStacks{i,1});
        end
        H.ZStack.c = bigStack;
        H.ZStack.reg = bigStack;
        
        %---------- Display First Image of Stack ----------%
        axes(H.axes1)
        H.image = imagesc(H.ZStack.c(:,:,1),[H.Min H.Max]); colormap(CubeHelix(256,0.5,-1.5,1.2,1.0))
        set(gca,'Box','on','XTick',[],'YTick',[])

        %---------- Send Info to GUI ----------%
        set(H.slider.frame,'Max',H.NFrames,'Min',1,'SliderStep',[1/H.NFrames 10/H.NFrames],'Value',1);
        set(H.slider.min,'Min',0,'Max',H.Max,'Value',H.Min);
        set(H.slider.max,'Min',0,'Max',H.Max,'Value',H.Max);
        set(H.rdio.zproj,'Value',0);
        set(H.rdio.roi,'Value',0);
        set(H.rdio.stack,'Value',1);
        H.ShowWhat = 'Stack';
        clc
        
    %---------- Filter ----------%
    function Filter(~,~)
    global H
    
        %----- Filter What -----%
        fw = get(H.popup.filtwhat,'String'); fw = fw{get(H.popup.filtwhat,'Value')};
        switch fw
            case 'Stack'
                ts = H.ZStack.c;
                set(H.rdio.zproj,'Value',0);
                set(H.rdio.roi,'Value',0);
                set(H.rdio.stack,'Value',1);
            case 'Proj'
                wp = get(H.popup.wproj,'String'); wp = wp{get(H.popup.wproj,'Value')};
                switch wp
                    case 'Mean'
                        ts = H.zproj.m;
                    case 'STD'
                        ts = H.zproj.s;
                    case 'CC'
                            if ~isfield(H.zproj,'c')
                                pp = gcp('nocreate');
                                if isempty(pp)
                                    parpool(H.npools); H.poolobj = gcp('nocreate');
                                end
                                H.zproj.c = CrossCorrImageKP(double(H.ZStack.c),H.npools,'CC','yes','no');
                        end
                        ts = H.zproj.c;
                end    
                set(H.rdio.zproj,'Value',1);
                set(H.rdio.roi,'Value',0);
                set(H.rdio.stack,'Value',0);
        end
                
        %----- Filter How -----%                             
        fh = get(H.popup.filthow,'String'); fh = fh{get(H.popup.filthow,'Value')};
        switch fh
            case 'Kalman'
                ts = Kalman_Stack_Filter(ts,.5);                                
            case 'Median'
                ts = medfilt3(ts,[3 3 3]);                
            case 'Wiener'
                parfor i = 1:size(ts,3)
                    ts(:,:,i) = wiener2(ts(:,:,i),5);             
                end          
            case 'Original'
                ts = H.ZStack.org;  
                H.zproj.m = mean(H.ZStack.c,3);
                H.zproj.s = std(double(H.ZStack.c),1,3);
                H.zproj = rmfield(H.zproj,'c');
                
        end        
        if size(ts,3) == 1
            ts = (ts - min(ts(:)))/(max(ts(:)) - min(ts(:)));
        end
        H.Max = max(ts(:));
        H.Min = min(ts(:));
        
         %----- Display First Image of Stack -----%
        axes(H.axes1)
        H.image = imagesc(ts(:,:,1),[H.Min H.Max]); colormap(CubeHelix(256,0.5,-1.5,1.2,1.0))
        set(gca,'Box','on','XTick',[],'YTick',[])

        %----- Send Info to GUI -----%
        if size(ts,3) > 1
            H.ZStack.c = ts;
        else
            switch wp
                case 'Mean'
                    H.zproj.m;
                case 'STD'
                    H.zproj.s = ts;
                case 'CC'
                    H.zproj.c = ts;
            end   
        end
        set(H.slider.frame,'Max',H.NFrames,'Min',1,'SliderStep',[1/H.NFrames 10/H.NFrames],'Value',1);
        set(H.slider.min,'Min',0,'Max',H.Max,'Value',H.Min);
        set(H.slider.max,'Min',0,'Max',H.Max,'Value',H.Max);
        
        H.ShowWhat = 'Stack';
        clc
        
        
%--------------- Find ROIs and Segmentation ---------------%
    
    %---------- Auto ROIs ----------%
    function AutoRoi(~,~)
    global H
        %----- Initiate -----%
        rt = get(H.popup.roihow,'String'); rt = rt{get(H.popup.roihow,'Value')};
        wproj = get(H.popup.roiwhat,'String'); wproj = wproj{get(H.popup.roiwhat,'Value')};
        p = get(H.edit.roiparam,'String');        
        switch wproj
            case 'Mean'
                Proj = H.zproj.m;
            case 'STD'
                Proj = H.zproj.s;
            case 'CC'
                if ~isfield(H.zproj,'c')
                    pp = gcp('nocreate');
                    if isempty(pp)
                        parpool(H.npools); H.poolobj = gcp('nocreate');
                    end
                    H.zproj.c = CrossCorrImageKP(double(H.ZStack.c),H.npools,'CC','yes','no');
                end
                Proj = H.zproj.c;                    
        end
        Proj = (Proj - min(Proj(:)))/(max(Proj(:)) - min(Proj(:)));
        
        %----- Get ROIs -----%
        switch rt
            case 'Laplace'
                if isempty(p)
                    clear p
                    p(1) = round2(mode(Proj(:)) + 2*std(Proj(:)),0.01);
                    p(2) = 100;
                    p(3) = 1000;
                    set(H.edit.roiparam,'String',num2str(p));
                else                
                    p = str2num(p);               
                end
                [H.roi.mask, laplace] = mthresh(Proj, p(1), p(2), p(3));

            case 'Threshold'         
                
                if isempty(p)
                    clear p
                    p(1) = round2(mode(Proj(:)) + 2*std(Proj(:)),0.01);
                    p(2) = 100;
                    p(3) = 1000;
                    set(H.edit.roiparam,'String',num2str(p));
                else                
                    p = str2num(p);               
                end                     
                
                I2 = imtophat(Proj,strel('disk',15));
                I3 = imadjust(I2);        
                mask2 = im2bw(I3,p(1));
                mask2 = imclearborder(mask2);
                H.roi.mask = bwareaopen(mask2, p(2),8);

            case 'ICA'

            case 'Correlation'
                %----- Threshold Image -----%                
                if isempty(p)Ima
                    clear p
                    p(1) = round2(mode(Proj(:)) + 2*std(Proj(:)),0.01);
                    p(2) = 100;
                    p(3) = 250;
                    p(4) = .15;
                    set(H.edit.roiparam,'String',num2str(p));
                elseif  length(str2num(p)) < 4      
                    clear p
                    p = str2num(p); 
                    p(4) = .15;
                    set(H.edit.roiparam,'String',num2str(p));
                else
                    p = str2num(p);               
                end
                I2 = imtophat(Proj,strel('disk',15));
                I3 = imadjust(I2);        
                mask2 = im2bw(I3,p(1));
                mask2 = imclearborder(mask2);
                H.roi.mask = bwareaopen(mask2, p(2),8);
                H.roi.reg = regionprops(H.roi.mask,'all');
                H.roi.con = bwconncomp(H.roi.mask, 8);
                [nx, ny, nt] = size(Proj);
                rankThresh = p(1);                
                minSize = p(2);
                maxSize = p(3);
                corrThresh = p(4);
                p(6) = 0;
                T=[];    
                H.roi.mask = zeros(nx,ny);
                nn = 0;
                for i = 1:H.roi.con.NumObjects 
                    set(H.edit.activity,'String',['Currently Dissecting ROI: ' num2str(i) ' of ' num2str(H.roi.con.NumObjects)]); pause(H.delay)                    
                    T(:,1) = H.roi.reg(i).PixelIdxList;
                    T(:,2) = Proj(T(:,1));
                    T = sortrows(T,2);
                    SI = T(:,1);                                                              
                    n = 0;  IDX = []; OSize = length(SI);
                    while ~isempty(SI)
                        n = n + 1; 

                        %--- Using Walking Algorythm ---%
                        [xi, yi] = ind2sub(size(Proj),SI(1));
                        if n == 1
                            BLOCKED = ones(nx,ny);
                            BLOCKED(H.roi.reg(i).PixelIdxList) = 0;                                             
                            [IDX, BLOCKED, newBLOCKED] = recWalkSubstitute(double(H.ZStack.c),BLOCKED, xi, yi, corrThresh, maxSize);
                        else 
                            corrThresh = p(4) - p(6)*n;
                            [IDX, BLOCKED, newBLOCKED] = recWalkSubstitute(double(H.ZStack.c), BLOCKED,xi, yi, corrThresh, maxSize);
                        end
                        if size(IDX,1) > minSize && size(IDX,1) < maxSize
                            nn = nn + 1;
                            for j = 1:size(IDX,1)
                                H.roi.mask(IDX(j,1), IDX(j,2) ) = nn;
                            end
                        end

                        %--- Get Rid of Blocked ---%
                        nB = find(newBLOCKED == 1);
                        for k = 1:length(nB)
                            J = find(nB(k) == SI);
                            if ~isempty(J) 
                                SI(J) = [];
                            end
                        end           

                    end
                    clear T SI n IDX OSize Blocked newBlocked
                end                                
                                 
        end 
        
        %----- Segment Mask -----%
%         H.roi,mask = pad(H.roi.mask,5);
        H.roi.map = H.roi.mask;
        H.roi.reg = regionprops(H.roi.mask,'all');
        H.roi.con = bwconncomp(H.roi.mask, 8);                

        %----- Display Results -----%        
        roitype = get(H.popup.roitype,'String'); roitype = roitype{get(H.popup.roitype,'Value')};
        axes(H.axes1)
        H.image = imagesc(Proj); hold on
        colormap(CubeHelix(256,0.5,-1.5,1.2,1.0))
        set(gca,'Box','on','XTick',[],'YTick',[]);
        for i = 1:H.roi.con.NumObjects            
            switch roitype
                case 'Scatter'                    
                    scatter(H.roi.reg(i).PixelList(:,1),H.roi.reg(i).PixelList(:,2),50,H.roiCM(i,:),'fill','MarkerEdgeColor','none');
                case 'T-Scatter'
                    scatter_patches(H.roi.reg(i).PixelList(:,1),H.roi.reg(i).PixelList(:,2),50,H.roiCM(i,:),'FaceAlpha',0.1,'EdgeColor','none');
                case 'Patch'                    
                    patch(H.roi.reg(i).ConvexHull(:,1),H.roi.reg(i).ConvexHull(:,2),'-','FaceColor',H.roiCM(i,:),'FaceAlpha',0.5,'EdgeColor','none');                    
                case 'Outline'
%                     contour(H.roi.reg(i).Image,1,'Color',H.roiCM(i,:));                    
                    plot(H.roi.reg(i).ConvexHull(:,1),H.roi.reg(i).ConvexHull(:,2),'-','Color',H.roiCM(i,:),'LineWidth',2);                    
            end            
        end 
        hold off
        
        %---------- Send Info to GUI ----------%
        H.Max = max(H.image(:)); H.Min = min(H.image(:));
        set(H.slider.frame,'Max',2,'Min',1,'SliderStep',[1 1],'Value',1);
        set(H.slider.min,'Min',0,'Max',H.Max,'Value',H.Min);
        set(H.slider.max,'Min',0,'Max',H.Max,'Value',H.Max);
        set(H.rdio.zproj,'Value',0);
        set(H.rdio.roi,'Value',1);
        set(H.rdio.stack,'Value',0);
        H.ShowWhat = 'ROI';
        set(H.edit.activity,'String','Finished finding ROIs');
        clc
        
    %---------- Manuel ROI ----------%
    function AddRoi(~,~)
    global H
    
        %----- Select Method -----%        
        H.roi.list = imellipse(H.axes1)
        
        %----- Save ROI -----%
        
        
        
        %----- Clean up GUI -----%
        set(H.rdio.addroi,'Value',0);        
    
    
    
        
    
    
        
        
    
    

            
            
            
            
            
            
