function [dat,ysize,xsize,nframes,metadata] = imreadND2(fname)

% FUNCTION [DAT,YSIZE,XSIZE,NFRAMES,METADATA] = IMREADND2(FNAME);
% 
% This function is written to load in an ND2 file taken using the Nikon
% Elements program.  It is specifically constructed in this case to import
% time stacks.  This function is based on the open microscopy core that was
% downloaded on 11/4/2014.
%
% DAT = Output data in an ysize x xsize x nframes array.
% *SIZE = Number of (x,y) pixels in movie
% NFRAMES = Number of imported frames in the movie
% METADATA = The recorded parameters imported from the Elements program
%
% Created by Joel Kralj.  2014-11-05
%
% Modified by Joel Kralj to enable multi-position movies.  If there is more
% than 1 xy position, dat will be a 4D array in the order, Y,X,t,FOV.

evalc('tmpdat = bfopen(fname);');
[nPos,~] = size(tmpdat);

for n = 1:nPos
    nframes(n) = length(tmpdat{n,1});
    tmp = tmpdat{n,1}(1);
    [ysize(n),xsize(n)] = size(tmp{1});

    metadata{n} = tmpdat{n,2};
    metadataKeys = metadata{n}.keySet().iterator();
    
    if length(tmpdat{1,1}) <= 2
        nframes = 1;
        dat = zeros(ysize,xsize,nframes,'uint16');
        dat(:,:,1,n) = double(tmp{1});
    else
        dat = zeros(ysize,xsize,nframes,'uint16');
        dat(:,:,:,n) = reshape(cell2mat(tmpdat{n,1}(1:nframes)),[ysize,xsize,nframes]);
    end
    
    
    indx = 0;
    for i=1:metadata{n}.size()
        try
            indx = indx+1;
            key{indx} = metadataKeys.nextElement();
            value{indx} = metadata{n}.get(key{indx});
        catch
        end
    end
    
    %% Sort values to more easily read in .txt file
    [key2,idx] = sort([key]);
    value2 = value(idx);
    
    txtFileName = [fname(1:end-4) '_metadata.txt'];
    fID = fopen(txtFileName,'w');
    for i = 1:indx
        fprintf(fID,'%s = %s\r\n', key2{i}, value2{i});
    end
    fclose(fID);
end