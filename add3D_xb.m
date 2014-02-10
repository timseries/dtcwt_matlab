function add3D_xb(filename,gain)

% function add3D(filename,gain)
% Add gain * the 3-D array of 32-bit floats from filename.3d
% to global array xb.

global xb

fid=fopen([filename '.3d'],'r');
if fid < 0,
   error(['Cannot open 3D file: ' filename '.3d for reading']);
end

sxb = fread(fid,3,'int32');

for f = 1:sxb(3),
   xb(:,:,f) = xb(:,:,f) + gain * fread(fid,[sxb(1) sxb(2)],'float32');
end

fclose(fid);
return;
