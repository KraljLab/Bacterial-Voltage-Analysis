function [illum,tForm] = find_illum(flist,primaryColor)
%
%  This function is designed to input a list of files all taken under
%  identical conditions, and then output a smoothed average of all the
%  data.  This is expecting files in a .nd2 format, though future versions
%  may incorporate handling .tif files.
%
%  illum = The illumination pattern - it will be of dimensions
%  (ysize,xsize,ncolors)
%  tForm = Transform to align the colors to the first color, if there is
%  more than one color
%
%  flist = a list of files generated from the matlab dir command.  This
%  function will use this variable to and use the names component of the
%  structure
%  primaryColor = if there is more than one color, this is used as the
%  stationary image to which all other colors are aligned
%
%  Written by Joel Kralj - 2018-05-02

rgAlign = 0;
nfiles = length(flist);
tForm = [];

for f = 1:nfiles
    fprintf(1,'calculating illumination for: %s\n',flist(f).name);
    evalc('nd2r = BioformatsImage(fullfile(pwd,flist(f).name));');
    
    if ~exist('primaryColor')
        tmp = strfind(nd2r.channelNames,'488');
        tmp = [tmp strfind(nd2r.channelNames,'Quad 2')];
        
        primaryColor = 2-mod(find(not(cellfun('isempty', tmp))),2);
    end
    StartFrame = round(nd2r.sizeT.*0.4);
    
    ncolors = nd2r.sizeC;
    ysize = nd2r.height;
    xsize = nd2r.width;
    
    tmpImg = zeros(ysize,xsize,10,ncolors);
    count = 0;
    for jj=StartFrame:StartFrame+9
        count = count + 1;
        for c = 1:ncolors
            tmpImg(:,:,count,c) = nd2r.getXYplane(c,1,jj); %Quad 1 i.e. tmp(:,:,1:2:end)
        end
    end
    
    if rgAlign == 0 && ncolors > 1
        toAlign = 1:ncolors;
        toAlign(primaryColor) = [];
            
        [optimizer,metric] = imregconfig('multimodal');
        stationary = adapthisteq(mat2gray(squeeze(mean(tmpImg(:,:,:,primaryColor),3))));
        for c = toAlign
            moving = adapthisteq(mat2gray(squeeze(mean(tmpImg(:,:,:,c),3))));
            tFormTmp = imregtform(moving,stationary,'similarity',optimizer,metric);
            tForm{c} = tFormTmp;
        end
        pixAdj = sum(abs(tForm{c}.T(:)))-3;
        if abs(pixAdj) < 30
            rgAlign = 1;
        else
            rgAlign = 0;
            tForm{c}.T = [1 0 0; 0 1 0; 0 0 1];
        end
    end
    
    % Start the matrix if this is the first time through
    if f == 1
        illumAll = zeros(ysize,xsize,ncolors,nfiles,'uint16');
    end
    
    % Capture the red and green average images
    for c = 1:ncolors
        size(imgaussfilt(mean(tmpImg(:,:,:,c),3),[45 45]));
        illumAll(:,:,c,f) = uint16(imgaussfilt(mean(tmpImg(:,:,:,c),3),[65 65]));
    end
end

% Broaden the illumination pattern by filtering
for c = 1:ncolors
    imgTmp = illumAll(:,:,c,:);
    illum(:,:,c) = imgaussfilt(imopen(mean(imgTmp,4),strel('ball',90,50)),[75 75]);
end

end
