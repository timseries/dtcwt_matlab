function spatial_transform
% SPATIAL_TRANSFORM 
% Performing 3-D interpolation on xa, using shift_x, shift_y and shift_z.
% The 3-D interpolation is carried out slice by slice, in order to avoid
% the huge amount of memory consumption.
global xa shift_x shift_y shift_z xa_before_reg

sxa = size(xa);
s2 = sxa(1); % no of rows
s1 = sxa(2); % no of columns
s3 = sxa(3); % no of frames

% Update shift_x shift_y and shift_z to the true shift instead of the shift
% wrt. the image, and also add the amount of matrix index to the motion
% field.

thx_2d = ones(s2,1) * [1:s1]; 
for k=1:s3, shift_x(:,:,k) = shift_x(:,:,k) * (s1/2) + thx_2d ; end
clear thx_2d
thy_2d = [1:s2]' * ones(1,s1); 
for k=1:s3, shift_y(:,:,k) = shift_y(:,:,k) * (s2/2) + thy_2d; end

for k=1:s3,
    shift_z(:,:,k) = shift_z(:,:,k)* (s3/2) + k ;
end

% Assigning zeros to the information being shifted from the outside of the dataset.
index = find(shift_y<1); shift_y(shift_y<1) = 1; 
index = union(index,find(shift_y>s2)); shift_y(shift_y>s2) = s2;

index = union(index,find(shift_x<1)); shift_x(shift_x<1) = 1;
index = union(index,find(shift_x>s1)); shift_x(shift_x>s1) = s1;

index = union(index,find(shift_z<1)); shift_z(shift_z<1) = 1;
index = union(index,find(shift_z>s3)); shift_z(shift_z>s3) = s3;

% To tackle the problem of huge memory consumption used by direct inter3,
% we divide the whole 3-D dataset into several blocks of data and carry out
% interp3 in these sub-blocks
for frame = 1:s3
    % For the current frame, find max shift in z direction
    max_shift_z = max(ceil(max(max(shift_z(:,:,frame)-frame))),0);
    % For the current frame, find min shift in z direction
    min_shift_z = min(floor(min(min(shift_z(:,:,frame)-frame))),0);
    z_region = max((frame+min_shift_z),1):min((frame+max_shift_z),s3);
   
    if (size(z_region,2)~=1)
        temp = interp3(xa_before_reg(:,:,z_region),shift_x(:,:,z_region),shift_y(:,:,z_region),shift_z(:,:,z_region)-z_region(1)+1);
        xa(:,:,frame) = temp(:,:,1+abs(min_shift_z)); 
    else 
        xa(:,:,frame) = interp2(xa_before_reg(:,:,frame),shift_x(:,:,frame),shift_y(:,:,frame));
    end
end

% xa = interp3(xa_before_reg,shift_x,shift_y,shift_z,'lin');
xa(index) = 0;
return