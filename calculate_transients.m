% This script will use the data extracted from extract_video_data.m and
% calculate calcium transients.  Adjust the directory variable as needed.
%
% Note this script is expecting an xls file to read out the conditions.
% In column 1 should be the movie, column 2 should be the gene knockout
% name, column3 should be the aminoglycoside concentration, and column 4
% should be the well name.

clear
close all

% Move to the directory of interest
cd('C:\VideoData');
load Long_term_corrected_extracted_data.mat

saveDir2 = [saveDir filesep 'IndividualPlots'];
if ~exist(saveDir2)
    mkdir(saveDir2)
end

%% Read in the experimental conditions from the xls file
xlsFile = 'Well_treatments.xlsx';
[~,~,rawData] = xlsread(xlsFile);

wellNames = rawData(1:end,4);

for f = 1:nfiles
    tmp1 = regexp(flist(f).name,'_','split');
    currWell = (tmp1{2}(1:2));
    if length(currWell) == 2
        currWell = [currWell(1) '0' currWell(2)];
    end
    currMov = find(contains(wellNames,currWell));
    data(f).movie = rawData{currMov,1};
    data(f).KO = rawData{currMov,2};
    data(f).gent = rawData{currMov,3};
    data(f).fname = flist(f).name;
    data(f).well = currWell;
end
clear rawData

%%
timeToPlot = 1:920; % Number of frames to plot
toPlot = 1;         % If you want to plot out images
for f = 1:nfiles
    % For each file, start with the intensity and masks
    intensG = intensGAll{f};
    [ncells, nframes] = size(intensG);
    
    tmp = timeAll{f};
    time = tmp(2,:);
    transT = find(time > KanAdd*3600,1);
    L = LGAll{f};
    
    % Define a title name for the figures
    titleName = [flist(f).name(end-5:end-4) ': ' data(f).KO ' - ' ...
        num2str(data(f).gent) ' ug/mL gent'];
    
    % Plot out a few random cells
    if toPlot
        try
            rIndx = randsample(ncells,5);
            figure
            plot(time(timeToPlot)/3600,intensG(rIndx,timeToPlot)')
            title([titleName ' : Random Traces'])
        catch
            disp('Could not find any cells...')
%             continue;
        end
    end
    
    % Plot out sanity check
    nToPlot = min(5,ncells);
    randIndx = randsample(ncells,nToPlot);
    
    if toPlot & nToPlot
        figure
        subplot(2,2,1:2)
        maskTmp = double(logical(L));
        maskTmp(1,1) = 2.5;
        imshowpair(medGAll{f},maskTmp)
        title(titleName)
        
        subplot(2,2,3:4)
        plot(time(timeToPlot)/3600,intensG(randIndx,timeToPlot)')
        line([time(transT(1)) time(transT(1))]/3600,[min(intensG(:)) .5*max(intensG(:))],'Color',[0 0 0])
        xlabel('Time (hrs)')
        ylabel('GCaMP6 fluorescence')
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_sanity_check.png'])
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_sanity_check.fig'])
    end
    
    % Attempt to quantify blinkiness of traces
    movMedTime = 51;
    movMedFrame = find(time(2,:)/60 > movMedTime,1);
    movStdTime = 31;
    movStdFrame = find(time(2,:)/60 > movStdTime,1);
    
    % Find the background, divide, and calculate moving standard deviation
    % for each cell
    bgTrace = movmedian(intensG,movMedFrame,2);
    intensG2 = bsxfun(@rdivide,intensG,bgTrace);
    intensG3 = movstd(intensG2',movStdFrame)';
    
    if toPlot
        try
        figure
        nplot = 9:12;
        figure
        subplot(2,1,1)
        plot(intensG(nplot,:)')
        subplot(2,1,2)
        plot(intensG3(nplot,:)')
        catch
        end
    end
    
    % Take the mean of the moving standard deviation across all cells
    meanStd = nanmean(intensG3,1);
    meanTot = nanmean(intensG,1);
    
    if toPlot
        figure
        plotyy(time(timeToPlot)/3600,meanTot(timeToPlot),...
            time(timeToPlot)/3600,meanStd(timeToPlot))
        title([titleName ' Mean and standard deviation of all traces'])
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2)  '_Mean_std_trace.png'])
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2)  '_Mean_std_trace.fig'])
    end
    
    % Calculate if there's a large transient
    blinkythresh = .1;
    clear peakAll
    for nActive = 1:ncells
        highTimes = find(intensG3(nActive,:) > blinkythresh);
        clear indx1 indx2 peakLength peak
        
        crossOvers = [find(diff(highTimes) > 1) length(highTimes)];
        nCrossOvers = length(crossOvers);
        
        for c = 1:nCrossOvers
            if c == 1
                indx1(c) = 1;
                indx2(c) = crossOvers(c) + 1;
            else
                indx1(c) = crossOvers(c-1)+1;
                indx2(c) = crossOvers(c) + 1;
            end
            peakLength(c) = indx2(c)-indx1(c);
        end
        
        lengthThresh = 40; % Length in frames of how long a peak needs to be
        goodPeaks = find(peakLength > lengthThresh);
        nGoodPeaks = length(goodPeaks);
        
        if nGoodPeaks > 0
            for c = 1:nGoodPeaks
                peak{c} = highTimes( indx1(goodPeaks(c)):indx2(goodPeaks(c))-1 );
            end
        else
            peak = 0;
        end
        
        peakAll{nActive} = peak;
        
    end
    
    %  Identify first time to transients
    didBlink = [];
    for c = 1:ncells
        didBlink(c)  = iscell(peakAll{c});
    end
    blinkers = find(didBlink);
    pctBlinkers = length(find(didBlink))/length(didBlink);
    
    indx = 0; firstBlink = [];
    for c = blinkers
        indx = indx+1;
        firstBlink(indx) = time(peakAll{1,c}{1,1}(1))/3600;
    end
    
    if toPlot & nToPlot
        figure
        set(gcf,'Position',[1281 701 360 290]);
        h = histogram(firstBlink,ceil(ncells/5),'Normalization','Probability');
        x = [time(transT)/3600 time(transT)/3600];
        y = [min(h.Values) max(h.Values)];
        line(x,y,'Color',[0 0 0])
        text(0.5,0.1,[num2str(pctBlinkers) '% blink'])
        title([titleName ' Start time of blinkiness'])
        xlabel('Blink initiation (hours)')
        ylabel('% blinking cells')
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_start_time_blinkys.png'])
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_start_time_blinkys.fig'])
    end
    
    % Look at time difference between AB addition and transient start
    try
        [yOut,xOut] = hist(firstBlink,round(ncells/3));
        fitBlinkG = fit(xOut',yOut','Gauss1');
        mostBlinky = fitBlinkG.b1;
        
        if toPlot
            figure
            histogram(firstBlink,round(ncells/3))
            hold on
            plot(fitBlinkG)
            hold off
            text(prctile(xOut,70),prctile(yOut,90),['\Delta = ' num2str(mostBlinky-time(transT)/3600,3) ' hrs'])
            saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_start_time_blinkys_with_fit.png'])
            saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_start_time_blinkys_with_fit.fig'])
        end
    catch
        disp(['Could not fit file ' num2str(f)]);
        mostBlinky = transT;
    end
    
    % Plot a few blinky vs non-blinky cells
    nToPlot = 9;
    nonBlinkers = find(~didBlink);
    
    randBli = randsample(length(blinkers),min(nToPlot,length(blinkers)));
    randNon = randsample(length(nonBlinkers),min(nToPlot,length(nonBlinkers)));
    
    meanStdBli = nanmean(intensG3(blinkers,:),1);
    meanTotBli = nanmean(intensG(blinkers,:),1);
    meanStdNon = nanmean(intensG3(nonBlinkers,:),1);
    meanTotNon = nanmean(intensG(nonBlinkers,:),1);
    
    if toPlot
        figure
        try
            subplot(2,3,1:2)
            plot(time(timeToPlot)/3600,intensG(blinkers(randBli),timeToPlot))
            subplot(2,3,3)
            plotyy(time(timeToPlot)/3600,meanTotBli(timeToPlot),...
                time(timeToPlot)/3600,meanStdBli(timeToPlot))
        catch
        end
        try
            subplot(2,3,4:5)
            if ~isempty(nonBlinkers)
                plot(time(timeToPlot)/3600,intensG(nonBlinkers(randNon),timeToPlot))
            end
        catch
        end
        subplot(2,3,6)
        plotyy(time(timeToPlot)/3600,meanTotNon(timeToPlot),...
            time(timeToPlot),meanStdNon(timeToPlot))
        title(titleName)
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_random_blinky_vs_nonblinky.png'])
        saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_random_blinky_vs_nonblinky.fig'])
    end
    
    growthTmp = growthGAll{f};
    growthTmpN = growthTmp./mean(growthTmp(1:10));
    growthTmpS = smooth(growthTmpN,21);
    growthRateS = smooth(diff(growthTmpS));

    if toPlot
        figure
        try
            h = histogram(firstBlink,round(ncells/5),'Normalization','Probability',...
                'EdgeAlpha',0.5, ...
                'FaceAlpha',0.5);
            ylabel('Cells blinking')
            yyaxis right
            plot(time(3:end)/3600,smooth(growthRateS(2:end),'sgolay'),...
                'Linewidth',1.5)
            ylabel('Growth rate')
            title(titleName)
            saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_blinking_vs_growth.png'])
            saveas(gca,[saveDir2 filesep num2str(f) '_' flist(f).name(1:2) '_blinking_vs_growth.fig'])
        end
    end
    
    growthRateAll{f} = growthRateS;
    meanStdAll{f} = meanStd;
    meanTotAll{f} = meanTot;
    pctBlinkersAll(f) = pctBlinkers;
    timeToStartBlinkAll(f) = mostBlinky-time(transT)/3600;
    xOutAll{f} = xOut;
    yOutAll{f} = yOut;
    
    clearvars -except *Dir* f flist *All nfiles toPlot KanAdd timeToPlot data ...
        KanRem
    close all
end

save Long_term_corrected_extracted_data.mat
       
        