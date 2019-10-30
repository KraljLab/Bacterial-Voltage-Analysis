function [L,ncells] = Hess6(img,SD,resizeScale,toPlot,minIntens)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Process the input arguments
if nargin == 1
    SD = 1.5;
    resizeScale = 2;
    toPlot = 1;
    minIntens = 150;
elseif nargin == 2
    resizeScale = 2;
    toPlot = 1;
    minIntens = 150;
elseif nargin == 3
    toPlot = 1;
    minIntens = 150;
elseif nargin == 4
    minIntens = 150;
end
%%
% Resize original image
i0 = imresize(double(img),resizeScale);

if resizeScale > 1
    i0 = medfilt2(i0,[3 3]);
end

% Run fibermetric to isolate curvature
i1 = fibermetric(i0,3.5*resizeScale,...
    'ObjectPolarity','bright',...
    'StructureSensitivity',SD);

i2 = im2bw(mat2gray(i1),graythresh(mat2gray(i1)));

% i3 = imopen(i2,strel('square',2*resizeScale));

i4 = bwlabel(xor(bwareaopen(i2,30*resizeScale^2), bwareaopen(i2,900*resizeScale^2)));
% i5 = imresize(i4,1/resizeScale,'nearest');
% Reduce a little bit
stats = regionprops(i4,i0,'MeanIntensity');
tooDim = find([stats.MeanIntensity] < minIntens);

%
Ltmp = i4;
for c = 1:length(tooDim)
    badPix = find(i4 == tooDim(c));
    Ltmp(badPix) = 0;
end
Ltmp = reshape(Ltmp,size(i0));

L = bwlabel(imresize(Ltmp,1/resizeScale,'nearest'));
ncells = max(L(:)); 

if toPlot
    mask = double(logical(L));
    mask(1,1) = 2.5;
    
    figure
    imshowpair(img,mask);
end

