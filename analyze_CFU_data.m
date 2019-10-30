clear
close all


%% Move to directory and load in .mat file with cell counts
activeDir = 'C:\CFUTest';
cd(activeDir)

load CFUCounts.mat

time = [0 1 2 4];   % Time points, in hours
kan = 30;           % Drug concentrations
genes = {'BW25113','atpA','atpB'};

%%
% Reshape the data to more easily compare conditoins and time points
CFU2 = reshape(CFU_counts',[4 3 3]);
[nTimes,nReps,nCons] = size(CFU2);

% Normalize the data to the initial counts and calculate means and
% standard deviations
%
% CFU(Time,Rep,Condition)
% CFUMean(Time,Condition)
for cc = 1:nCons
    CFUNorm(:,:,cc) = CFU2(:,:,cc)./repmat(squeeze(CFU2(1,:,cc)),[4 1]);
    CFUNormMean(:,cc) = squeeze(mean(CFUNorm(:,:,cc),2));
    CFUNormStd(:,cc) = squeeze(std(CFUNorm(:,:,cc),[],2));
    CFUMean(:,cc) = squeeze(mean(CFU2(:,:,cc),2));
    CFUStd(:,cc) = squeeze(std(CFU2(:,:,cc),[],2));
end

% Plot out CFUs with error bars on a log scale
figure
cOrder = get(gca,'ColorOrder');
hold on
for cc = 1:nReps
    hh(cc) = semilogy(time,CFUNormMean(:,cc),'Color',cOrder(cc,:));
    errorbar(time,CFUNormMean(:,cc),CFUNormStd(:,cc)/sqrt(3),'Color',cOrder(cc,:))
end
set(gca,'YScale','log','YLim',[5e-6 1])
legend(hh,genes,'Location','SouthWest')
title([num2str(kan) ' \mug/mL kanamycin'])


saveas(gca,'FigureName.png')
saveas(gca,'FigureName.fig')

%% Test signficance at different conditions
timePt1 = 1; timePt2 = 4;
cond1 = 1;
cond2 = 1;

test1 = squeeze(CFU2(timePt1,:,cond1));
test2 = squeeze(CFU2(timePt2,:,cond2));

[h1,p1] = ttest2(test1,test2,'vartype','unequal')

