% This script is to extract single cell traces from E. coli taken on a
% Nikon Ti2 widefield microscope.  It assumes that you have 2 channels, red
% and green.

clear
close all

% Point to the active directory and define code conditions
activeDir = 'C:\VideoData';
cd(activeDir)

% Find files in the directory
flist = dir('*.nd2');
nfiles = length(flist);

% Define characteristics of scritp
toPlot = 0;         % Print out uncorrected movie?
writeAVI = 1;       % Write a truncated AVI for easing viewing
checkFocus = 0;     % Will search for frames that lose focus 
framesToKeep = 0;   % Set to 0 to go to maximum # frames
subBack = 1;        % Subtract background?
corrIllum = 1;      % Correct illumination?
KanAdd = 0.08;      % Time of kan addition in hours
KanRem = 0;         % Time of kan removal (in hours)
alignRG = 1;        % Align the red and green channels


%%
saveDir = [activeDir filesep 'Results_' num2str(date)];
codeDir = [saveDir filesep 'mfiles'];
if ~exist(saveDir)
    mkdir(saveDir)
    mkdir(codeDir)
end

% Make a copy of this code for posterity
try
    copyfile('*.m',codeDir);
catch
    disp('No mfiles in directory')
end


%% Find uneven illumination
if ~exist([saveDir filesep 'illum.mat'],'file')
    disp('Finding illumination of red and green')
    
    [illumTmp, tFormRG] = find_illum(flist); 
    illumG = uint16(illumTmp(:,:,2));
    illumR = uint16(illumTmp(:,:,1));
    
    save([saveDir filesep 'illum.mat'],'illumG','illumR','tFormRG');
    
    clearvars -except illumG illumR *Dir Kan* flist nfiles toPlot writeAVI...
        checkFocus framesToKeep subBack corrIllum tFormRG alignRG
else
    load([saveDir filesep 'illum.mat'])
end

figure
imshowpair(illumG,illumR)


%% Load in data
% Setup red to green alignment
[optimizer,metric] = imregconfig('multimodal');
[~,Reg] = imregister(illumR,illumG,'similarity',optimizer,metric);

% Loop through files
for f = 1:nfiles
    % Load in the video data
    disp([num2str(f) ' of ' num2str(nfiles) ': ' flist(f).name]);
    tic
    evalc('dat = imreadND2(flist(f).name);');
    toc
    
    %Extract the timestamp
    [time, Fs, camera] = extract_ND2_timestamp([flist(f).name(1:end-4) '_metadata.txt']);
    
    % If you want to truncate time points, this is where it happens
    if framesToKeep
        try
            time(:,framesToKeep+1:end) = [];
            R = dat(:,:,1:2:2*framesToKeep);
            G = dat(:,:,2:2:2*framesToKeep);
        catch
            nMax = 2*floor(size(dat,3)/2);
            time(:,nMax+1:end) = [];
            R = dat(:,:,1:2:2*nMax);
            G = dat(:,:,2:2:2*nMax);
        end
    else
        R = dat(:,:,1:2:end);
        G = dat(:,:,2:2:end);
    end
    [ysize,xsize,nFrames] = size(G);
    clear dat
    
    % Correct the illumination for each frame
    if corrIllum
        disp('   Dividing by uneven illumination...');
        GI = zeros(size(G),'uint16');
        RI = zeros(size(R),'uint16');
        for t = 1:nFrames
            tmp = mean(illumG(:))*double(G(:,:,t))./double(illumG);
            GI(:,:,t) = uint16(tmp);
            tmp = mean(illumR(:))*double(R(:,:,t))./double(illumR);
            RI(:,:,t) = uint16(tmp);
        end
        
        figure
        imshowpair(mean(GI,3),mean(RI,3))
        title('Illumination corrected mean images')
    else
        GI = G;
        RI = R;
    end
    clear G R
    
    
    % See if the focus drifts, this is useful if 1 or a few frames lose
    % focus during the movie
    if checkFocus
        disp('   Checking for out of focus frames...');
        for j = 1:nFrames
            tmp = double(GI(:,:,j));
            Gstd(j) = std(tmp(:));
        end
        dGstd = diff(Gstd);
        diffStd = std(dGstd(find(dGstd > prctile(dGstd,40)) & find(dGstd < prctile(dGstd,60))));
        
        badFrames = [];
        transDown = find(dGstd < mean(dGstd)-3*diffStd);
        for t = 1:length(transDown)-1
            try
                transUp(t) = find(dGstd(transDown(t):end) > mean(dGstd) + 3*diffStd,1)...
                    +transDown(t) + 1;
            catch
                try
                    transUp(t) = min(transDown(t+1),nFrames);
                catch
                    transUp(t) = nFrames;
                end
            end
            badFrames = unique([badFrames transDown(t):transUp(t)]);
        end
        disp(['   Removing ' num2str(length(badFrames)) ' bad frames...']);
        GI(:,:,badFrames) = [];
        RI(:,:,badFrames) = [];
        time(:,badFrames) = [];
        [~,~,nFrames] = size(GI);
    end
    
    
    % Display the entire movie
    cMaxG = prctile(GI(:),99.95);
    cMinG = (min(GI(:)) + .02*cMaxG);
    cMaxR = prctile(RI(:),99.95);
    cMinR = (min(RI(:)) + .02*cMaxR);
    
    if toPlot
        figure
        for t = 1:5:nFrames
            evalc('imshowpair(GI(:,:,t),RI(:,:,t))');
            drawnow
        end
    end
    
    % Align the red movies to the green
    if alignRG
        for t = 1:nFrames
            RI(:,:,t) = imwarp(RI(:,:,t),tFormRG{1},'Outputview',Reg);
        end
    end
    
    
    % Correct for drift and show images
    disp('   Correcting drift...')
    corrected = 0;
    yNew = 1:size(GI,1);
    xNew = 1:size(GI,2); 
    G2 = GI; R2 = RI;
    cIndx = 0;
    while ~corrected
        cIndx = cIndx+1;
        if cIndx <= 3
            trunc = 4;
        else
            trunc = 2;
        end
        tic
        [G2,R2,xShift,yShift] = correct_jitter(G2,R2,10,trunc);
        xNew = 1 + ceil(max(xShift)):xsize + ceil(min(xShift));
        yNew = 1 + ceil(max(yShift)):ysize + ceil(min(yShift));
               
        G2 = G2(yNew,xNew,:);
        R2 = R2(yNew,xNew,:);
        xsize = length(xNew);
        ysize = length(yNew);
        toc
        
        pixAdjust = sum(abs(diff(xShift))) + sum(abs(diff(yShift)));
        disp(['Curr shift err: ' num2str(pixAdjust)])
        if pixAdjust < nFrames*.025 || cIndx == 5
            corrected = 1;
        end
    end
    
    clear G R
    
    figure
    evalc('imshowpair(mean(G2,3),mean(R2,3));');

    
    % Try subtracting background
    if subBack
        for k = 1:nFrames
            R3(:,:,k) = R2(:,:,k) - imgaussfilt(imopen(R2(:,:,k),strel('disk',30)),[25 25]);
            G3(:,:,k) = G2(:,:,k) - imgaussfilt(imopen(G2(:,:,k),strel('disk',30)),[25 25]);
        end
        G3 = G3 - min(G3(:)) + 50;
        R3 = R3 - min(R3(:)) + 50;
    else
        G3 = G2;
        R3 = R2;
    end
    
    clear G2 R2
    
    
    % Write AVI if you asked for it
    if writeAVI
        colorImg = zeros(ysize*3,xsize,3);
        v = VideoWriter([flist(f).name(1:end-4) '.avi'],'Motion JPEG AVI'); %#ok<TNMLP>
        open(v);
        if nFrames < 300
            sp = 1;
        elseif nFrames < 800
            sp = 3;
        else
            sp = 5;
        end
        figure
        set(gcf,'Position',[680 165 1000 800])
        for t = 1:sp:nFrames
            rTmp = (double(R3(:,:,t)) - double(cMinR))/double(cMaxR);
            gTmp = (double(G3(:,:,t)) - double(cMinG))/double(cMaxG);
            % Green image
            colorImg(1:ysize,:,:) = repmat(gTmp,[1 1 3]);
            colorImg(ysize+1:2*ysize,:,:) = repmat(rTmp,[1 1 3]);
            colorImg(2*ysize+1:end,:,1) = rTmp;
            colorImg(2*ysize+1:end,:,2) = gTmp;
            
            currTime = datestr(seconds(time(2,t)),'HH:MM');
            
            if time(2,t) > 3600*24
                currTime(1:2) = num2str(24+str2num(currTime(1:2)));
            end
            evalc('imshow(colorImg)');
            text(10,10,[currTime],'Color',[1 1 1])
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        close(v);
    end
    
    
    % Look for cell growth by intensity
    if strcmpi(camera,'Andor')
        growThresh = 2000;
        brightThresh = 1800;
    else
        growThresh = 300;
        brightThresh = 200;
    end
    
    % Calculate percent of bright pixels (ie, with cells)
    nPix = length(yNew)*length(xNew);
    for t = 1:nFrames
        pctBrightG(t) = length(find(G3(:,:,t) > growThresh))/nPix;
        pctBrightR(t) = length(find(R3(:,:,t) > growThresh))/nPix;
    end
    
    % Segment and extract traces
    startFrame = min(601,nFrames-20);
    LG = Hess6(mean(G3(:,:,startFrame:1:startFrame+10),3),1.5,2,1,brightThresh);
    LR = Hess6(mean(R3(:,:,startFrame:1:startFrame+10),3),1.5,2,1,brightThresh);
    ncellsG = max(LG(:));
    ncellsR = max(LR(:));
    
    intensG = extract_traces(LG,G3(:,:,:),ncellsG);
    intensR = extract_traces(LG,R3(:,:,:),ncellsG);
    
    % Plot out a few random traces as a sanity check
    L2 = double(logical(LG));
    L2(1,1) = 2.5;
    nToPlot = min(5,ncellsG);
    randIndx = randsample(ncellsG,min(ncellsG,5));
    
    colorimg = zeros(ysize,xsize,3);
    colorimg(:,:,1) = 1.7*mat2gray(mean(R3,3));
    colorimg(:,:,2) = 1.7*mat2gray(mean(G3,3));
    colorimg(:,:,3) = 1.3*mat2gray(L2);
    
    figure
    subplot(2,2,1:2)
    evalc('imshow(colorimg)');
    title(flist(f).name,'interpreter','none')
    
    try
        subplot(2,2,3)
        plot(time(2,:)/3600,intensG(randIndx,:)')
        subplot(2,2,4)
        plot(time(2,:)/3600,intensR(randIndx,:)')
        if KanAdd
            subplot(2,2,3)
            yLims = get(gca,'YLim');
            line([KanAdd KanAdd],[yLims(1) yLims(2)],'color',[1 0 0])
            subplot(2,2,4)
            yLims = get(gca,'YLim');
            line([KanAdd KanAdd],[yLims(1) yLims(2)],'color',[1 0 0])
        end
        if KanRem
            subplot(2,2,3)
            yLims = get(gca,'YLim');
            line([KanRem KanRem],[yLims(1) yLims(2)],'color',[0 1 0])
            subplot(2,2,4)
            yLims = get(gca,'YLim');
            line([KanRem KanRem],[yLims(1) yLims(2)],'color',[0 1 0])
        end
    catch
        disp('Could not plot data...');
    end
    xlabel('Time (hrs)')
    ylabel('GCaMP Intensity')
    saveas(gca,[saveDir filesep num2str(f) '_sanity_check.png'])
   
    
    % Save out important variables
    intensGAll{f} = intensG;
    intensRAll{f} = intensR;
    LGAll{f} = LG;
    LRAll{f} = LR;
    timeAll{f} = time;
    medGAll{f} = mean(G3(:,:,round(nFrames/2):5:round(nFrames/2)+100),3);
    maskGAll{f} = mean(G3(:,:,startFrame:1:startFrame+10),3);
    medRAll{f} = mean(R3(:,:,round(nFrames/2):5:round(nFrames/2)+100),3);
    maskRAll{f} = mean(R3(:,:,startFrame:1:startFrame+10),3);
    growthGAll{f} = pctBrightG;
    growthRAll{f} = pctBrightR;
    
    clearvars -except *All activeDir saveDir flist nfiles f toPlot KanAdd...
        KanRem writeAVI checkFocus framesToKeep subBack corrIllum illum* ...
        tFormRG alignRG optimizer metric Reg
    close all
    
    
    % Save out variables to HDD
    save Long_term_corrected_extracted_data
end

