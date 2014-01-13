function savetif(vol, filename)
% SAVETIF saves volumetric data as a stacked tiff image

imwrite(vol(:,:,1), filename, 'WriteMode', 'overwrite');
for i = 2:size(vol,3);
    imwrite(vol(:,:,i), filename,'WriteMode','append');
end
