function [v,info] = readtif(tif)
% Converts a 3D stacked *.tiff file into a MATLAB mstrix
%
% SYNTAX:
%   v = readtif(tif);
%

info = imfinfo(tif,'tiff');
v = zeros(info(1).Width, info(1).Height, length(info)); 
for i = 1:length(info);
    v(:,:,i) = imread(tif, 'index', i);
end
