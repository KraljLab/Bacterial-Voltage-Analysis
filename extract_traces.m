function cellIntens = extract_traces(L,dat,numcells)

% FUNCTION CELLINTENS = EXTRACT_TRACES(L, DAT, NUMCELLS)
%
% This function is designed to apply a mask (L) generated by the
% segment_gcamp function and apply it to a movie (dat) to extract the
% time trace of each cell.
%
% CELLINTENS = The extracted traces
% L = Mask generated by segment_gcamp.m
% DAT = Movie data
% NUMCELLS = Number of cells output by segment_gcamp.m
% 
% Written by Joel Kralj. 2014-12-15

[~,~,nframes] = size(dat);

cellIntens = zeros(numcells,nframes);
for j = 1:nframes
    cellIntens(:,j) = cell2mat(struct2cell(regionprops(L,dat(:,:,j),'MeanIntensity')));
end
% cellIntens(1,:) = [];