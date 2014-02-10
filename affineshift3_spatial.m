function affineshift3_spatial( avec, avlevel )
% function affineshift3(avec, avlevel, mode) 
% Create 3-D shift vector (N.B shift vector contains 3 element which represent the shift along each of the 3 directions) from the affine vectors in avec(:,:,:,12)
% It is assumed that xi goes from -1 to +1 across the full width of the
% image (i.e. to the outer edges of the pixels).
% avec is the affine parameters obtained from the motion estimation, and
% avlevel is the level of the avec.

global shift_x shift_y shift_z
global band_size

% Increase size of avec by 1 if dimension is odd.
sav = size(avec);
if rem(sav(1),2) == 1
    t = 1:sav(1);
    avec = 0.5 * (avec([1 t],:,:,:) + avec([t sav(1)],:,:,:));
end
if rem(sav(2),2) == 1
    t = 1:sav(2);
    avec = 0.5 * (avec(:,[1 t],:,:) + avec(:,[t sav(2)],:,:));
end
if rem(sav(3),2) == 1
    t = 1:sav(3);
    avec = 0.5 * (avec(:,:,[1 t],:) + avec(:,:,[t sav(3)],:));
end

% Upsample the affine vector in two stages (with the purpose of reducing
% calculation time and memory usage)
% avec will be upsampled to the scale used in lev2, then generate motion
% field, finally upsample the lev2 motion field to acquire full resolution
% motion filed.
avec_upsample_to_lev = 2; 
  
% s1, s2 and s3 contain size of the scale that avec will be upsampled to
s1 = band_size(avec_upsample_to_lev,1)/2;
s2 = band_size(avec_upsample_to_lev,2)/2;
s3 = band_size(avec_upsample_to_lev,3)/2;

% Initialise the shift vectors
shift_y = zeros(s1,s2,s3);
shift_x = zeros(s1,s2,s3);
shift_z = zeros(s1,s2,s3);

% row / column / frame index distances are used to form the coordinate
% indexing when we generate shift vectors of affine parameters
row_index_distance = 2/band_size(1,1)*(2^avec_upsample_to_lev); % N.B band_size(1,:) gives original image size
column_index_distance = 2/band_size(1,2)*(2^avec_upsample_to_lev);
frame_index_distance = 2/band_size(1,3)*(2^avec_upsample_to_lev);

av = avec(:,:,:,1);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_y = shift_y + av;

av = avec(:,:,:,2);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_x = shift_x + av;

av = avec(:,:,:,3);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_z = shift_z + av;

%mat_index_2d = ((2*[1:s1]'-s1-1) / s1) * ones(1,s2); %#ok<NBRAK>
mat_index_2d = [-(s1-1)/2*row_index_distance:row_index_distance:(s1-1)/2*row_index_distance]' * ones(1,s2); %#ok<NBRAK>
mat_index = zeros(s1,s2,s3); %initialise mat_index
for k = 1:s3,
    mat_index(:,:,k) = mat_index_2d;
end % Now mat_index is actually yi

av = avec(:,:,:,4);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_y = shift_y + av.*mat_index;

av = avec(:,:,:,5);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_x = shift_x + av.*mat_index;

av = avec(:,:,:,6);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_z = shift_z + av.*mat_index;

%mat_index_2d = ones(s1,1)*((2*[1:s2]-s2-1) / s2); %#ok<NBRAK>
mat_index_2d = ones(s1,1) * [-(s2-1)/2*column_index_distance:column_index_distance:(s2-1)/2*column_index_distance];
for k = 1:s3,
    mat_index(:,:,k) = mat_index_2d;
end % Now mat_index is actually xi
av = avec(:,:,:,7);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_y = shift_y + av.*mat_index;
    
av = avec(:,:,:,8);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_x = shift_x + av.*mat_index;
    
av = avec(:,:,:,9);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_z = shift_z + av.*mat_index;

for k = 1:s3,
    %mat_index(1:s1,1:s2,k) = (2*k-s3-1) / s3;
    mat_index(1:s1,1:s2,k) = -(s3-1)/2*frame_index_distance + (k-1)*frame_index_distance;
end
av = avec(:,:,:,10);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_y = shift_y + av.*mat_index;

av = avec(:,:,:,11);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_x = shift_x + av.*mat_index;

av = avec(:,:,:,12);
for i = 1:(avlevel-avec_upsample_to_lev)
    av = upsample3(av,2);
    % Check if the coarser level subband is exactly twice as big as the
    % finer level subband. If not, discard 2 outer most coefs.
    if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
        av = av(2:end-1,:,:,:);
    end
    if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
        av = av(:,2:end-1,:,:);
    end
    if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
        av = av(:,:,2:end-1,:);
    end
end
shift_z = shift_z + av.*mat_index;

% Now upsample the motion field to full resolution, i.e. shift vectors will
% be upsampled to the original image resolution.
shift_y = upsample3(shift_y,4); % Now the whole shif_y is obtained
shift_x = upsample3(shift_x,4); % Now the whole shif_x is obtained
shift_z = upsample3(shift_z,4); % Now the whole shif_z is obtained
if band_size(1,1)/2 ~= band_size(2,1)
    shift_y = shift_y(3:end-2,:,:);
    shift_x = shift_x(3:end-2,:,:);
    shift_z = shift_z(3:end-2,:,:);
end
if band_size(1,2)/2 ~= band_size(2,2)
    shift_y = shift_y(:,3:end-2,:);
    shift_x = shift_x(:,3:end-2,:);
    shift_z = shift_z(:,3:end-2,:);
end
if band_size(1,3)/2 ~= band_size(2,3)
    shift_y = shift_y(:,:,3:end-2);
    shift_x = shift_x(:,:,3:end-2);
    shift_z = shift_z(:,:,3:end-2);
end

return;