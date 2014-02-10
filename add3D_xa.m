function add3D_xa(filename,gain)

% function add3D(filename,gain)
% Add gain * the 3-D array of 32-bit floats from filename.3d
% to global array xa.

global xa

fid=fopen([filename '.3d'],'r');
if fid < 0,
   error(['Cannot open 3D file: ' filename '.3d for reading']);
end

sxa = fread(fid,3,'int32');

for f = 1:sxa(3),
   xa(:,:,f) = xa(:,:,f) + gain * fread(fid,[sxa(1) sxa(2)],'float32');
end

fclose(fid);
return;
