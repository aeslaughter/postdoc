function saveinr(vol, filename)
% SAVEINR saves volumetric data as an *.inr file

% Write the inr file (iso2mesh)
fid = fopen(filename,'wb');
btype ='unsigned fixed';
dtype ='uint8';
bitlen = 8;

% Write the header
header=sprintf(['#INRIMAGE-4#{\nXDIM=%d\nYDIM=%d\nZDIM=%d\nVDIM=1\nTYPE=%s\n' ...
  'PIXSIZE=%d bits\nCPU=decm\nVX=1\nVY=1\nVZ=1\n'],size(vol),btype,bitlen);
header=[header char(10*ones(1,256-4-length(header))) '##}' char(10)];
fwrite(fid,header,'char');

% Write the data
fwrite(fid,vol,dtype);
fclose(fid);