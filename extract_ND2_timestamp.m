function [tOut, Fs, camera, nchannels, Laser] = extract_ND2_timestamp(txtName)
%
% This function is designed to open a text file generated reading in an ND2
% file with imreadND2_write_metadata and calculate out the corresponding
% timing of each frame.
%
% tOut is a 2 column matrix with the frame number in the first column, and
% the time stamp in the second.  Fs is the calculated sampling frequency.


%% Load in the metadata file
text = fileread(txtName);

%% Search for timestamps
p1 = '(\d+)';
p2 = '(\d)';
p3 = '(\d+)';
p4 = '((+|-)\d+)';
expr = ['timestamp #' p1 ' = ' p2 '.' p3 'e' p4];
t = regexp(text,expr,'tokens');

    
for q = 1:length(t)
    tmpVar = cellfun(@str2num,t{1,q});
    frNum(q) = tmpVar(1);
    nDiv = length(t{1,q}{3});
    tstamp(q) = (tmpVar(2) + (tmpVar(3)/10^(nDiv)))*10^(tmpVar(4));
end

%% Calculate out timing
[X,I] = sort(frNum);
Y = tstamp(I);

tOut = [X; Y];
Fs = 1/mean(diff(Y));

%% Find how many channels
try
    p5 = '(\d)';
    expr2 = ['Picture Planes = ' p5];
    t2 = regexp(text,expr2,'tokens');
    nchannels = cellfun(@str2num,t2{1});
catch
    nchannels = 1;
end

if strfind(text,'Andor')
    camera = 'Andor';
elseif strfind(text,'Flash')
    camera = 'Flash4';
else
    camera = 'Unknown';
end

%% Find out laser intensities
% p6 = '(\d+)'; % Laser color
% p7 = '(\d)'; % Exposure number
% p8 = '(\d+)'; % Intensity pct
% p9 = '(\d)';
% p10 = '(On|Off)'; % On or off
% expr3 = ['ExW:' p6 '; Power #' p7 ' = ' p8 '.' p9 '; ' p10];
% 
% try
%     t3 = regexp(text,expr3,'tokens');
%     indx = 0;
%     for q = 1:length(t3)
%         if ~isempty(find(ismember(t3{q}, 'On')))
%             indx = indx+1;
%             LasersOn(indx) = q;
%         end
%     end
%     
%     field1 = 'wavelength'; val1 = cell(nchannels,1);
%     field2 = 'intensity'; val2 = cell(nchannels,1);
%     Laser = struct(field1,val1,field2,val2);
%     for q = 1:length(LasersOn)
%         ExpNum = cellfun(@str2num,(t3{LasersOn(q)}(2)));
%         
%         Laser(ExpNum).intensity = cellfun(@str2num,t3{LasersOn(q)}(3)) + .1*cellfun(@str2num,t3{LasersOn(q)}(4));
%         Laser(ExpNum).wavelength = cellfun(@str2num,(t3{LasersOn(q)}(1)));
%     end
%     
% catch
%     t3 = regexp(text,['ExW:' p6 '; Power = ' p8 '.' p9 '; ' p10],'tokens');
% indx = 0;
%     for q = 1:length(t3)
%         if ~isempty(find(ismember(t3{q}, 'On')))
%             indx = indx+1;
%             LasersOn(indx) = q;
%         end
%     end
%     
%     field1 = 'wavelength'; val1 = cell(nchannels,1);
%     field2 = 'intensity'; val2 = cell(nchannels,1);
%     Laser = struct(field1,val1,field2,val2);
%     Laser(1).intensity = cellfun(@str2num,t3{LasersOn(1)}(2)) + .1*cellfun(@str2num,t3{LasersOn(1)}(3));
%     Laser(1).wavelength = cellfun(@str2num,(t3{LasersOn(1)}(1)));
    
Laser = 0;
% end
end