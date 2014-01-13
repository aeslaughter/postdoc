function v = tif2inr(tif, varargin)
% Converts a 3D stacked *.tiff file into a *.inr file
%
% This file was generated using the iso2mesh library source code.
%
% SYNTAX:
%   tif2inr(tif);
%   tif2inr(tif,inr)
%   v = tif2inr(...);
%
% DESCRIPTION:
%   tif2inr(tif) converts the *.tif file (tif) into a similarily named
%       *.inr file
%   tif2inr(tif,inr) converts the *.tif file into the specified *.inr file
%   v = tif2inr(...) same as above but returns the volume data

% Set the filename
if nargin == 1;
    [p,f,~] = fileparts(tif);
    inr = fullfile(p, [f, '.inr']);
else
    inr = varargin{1};
end

% Read the tif file into a volume array
[v, info] = readtif(tif);

% Write the inr file (iso2mesh)
fid = fopen(inr,'wb');
btype ='unsigned fixed';
dtype ='uint8';
bitlen = 8;

% Write the header
header=sprintf(['#INRIMAGE-4#{\nXDIM=%d\nYDIM=%d\nZDIM=%d\nVDIM=1\nTYPE=%s\n' ...
  'PIXSIZE=%d bits\nCPU=decm\nVX=1\nVY=1\nVZ=1\n'],size(v),btype,bitlen);
header=[header char(10*ones(1,256-4-length(header))) '##}' char(10)];
fwrite(fid,header,'char');

% Write the data
fwrite(fid,v,dtype);
fclose(fid);