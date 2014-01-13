function data = randvol(varargin)
%RANDVOL Creates a smooth random porous volume
%
% SYNTAX:
%   data = randvol;
%   data = randvol(n);
%
% DESCRIPTION:
%   data = randvol creates a 10 pixel cube random volume
%   data = randvol(n) creates an n pixel cube random volume
%   data = randvol(x,y,z) creates a volume with specified dimensions

% Determine the cube dimensions
if nargin == 0;
    x = 10; y = 10; z = 10;
elseif nargin == 1;
    N = varargin{1};
    x = N; y = N; z = N;
elseif nargin == 3;
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
else
    error('Input not reconginzed!');
end

% Create the random data
data = rand(x,y,z);

% Smooth the data
data = smooth3(data,'box',5);

% Set the data to 0 or 1
data(data<=0.5) = 0;
data(data>0.5) = 1;

% Convert data to unsigned 8-bit integer
data = uint8(data);

