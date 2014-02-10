function save3D(filename)

% function save3D(filename)
% Save the global 3-D array x in filename.3d using 32-bit floats.

global x

sx = size(x);
if length(sx) ~= 3,
   error('Not a 3-D array');
end

fid=fopen([filename '.3d'],'w');
if fid < 0,
   error(['Cannot open 3D file: ' filename '.3d for writing']);
end

disp(['Saving 3D data to ' filename '.3d'])
fwrite(fid,sx,'int32');
fwrite(fid,x,'float32');

fclose(fid);
return;



