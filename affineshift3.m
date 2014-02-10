function [shift_x shift_y shift_z] = affineshift3( avec, lev, avlevel )
% function [shift_x shift_y shift_z] = affineshift3(avec, lev, avlevel)
% The input parameters of the function: 
% avec => the affine parameters
% lev  => the level at which the shift vectors are to be generated
% avlevel => the level of the avec
% Create 3-D shift vector (N.B shift vector contains 3 element which represent the shift along each of the 3 directions) from the affine vectors in avec(:,:,:,12)
% It is assumed that xi goes from -1 to +1 across the full width of the
% image (i.e. to the outer edges of the pixels). 
% E.g. | x  x  x  x | x  x  x  x |      LEV n
%      |  o      o  |  o      o  |      LEV n+1
%     -1            0            +1
%
% avec is the set of affine parameters at current avlevel. We aim to
% upsample or downsample our avec from scale of avlevel to scale of lev. Then, using the
% up / down sampled affine vectors, generate the shift vectors at current
% lev.

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


av = avec;
sav = size(av);

if avlevel > lev
    av = zeros([band_size(lev,:)/2,12]); % initialise av
    for k = 1:sav(4) %note sav(4) should equal to 12, indicating 12 affine parameters
        av_temp = avec(:,:,:,k);
        for i = 1:(avlevel-lev)
         %   av = zeros([band_size(avlevel-i+1,:),12]); % initialise av
            % Upsample the affine vectors avec
            av_temp = upsample3(av_temp,2);
            % Check if the coarser level subband is exactly twice as big as
            % the finer level subband. If not, discard the 2 outer most
            % coefs.
            if band_size(avlevel-i,1)/4 ~= band_size(avlevel-i+1,1)/2
                av_temp = av_temp(2:end-1,:,:);
            end
            if band_size(avlevel-i,2)/4 ~= band_size(avlevel-i+1,2)/2
                av_temp = av_temp(:,2:end-1,:);
            end
            if band_size(avlevel-i,3)/4 ~= band_size(avlevel-i+1,3)/2
                av_temp = av_temp(:,:,2:end-1);
            end
        end
        av(:,:,:,k) = av_temp;
    end
elseif avlevel < lev     
    % If the coarser subband size is not a multiple of the finer
    % subband size, then we append the coefs by repeating edges
    % before the downsampling.
    
    for i = 1:(lev-avlevel)
        % check if the coarser level subband is exactly twice as big as the
        % finer level subband. If not, do appending.
        if band_size(avlevel+i-1,1)/4 ~= band_size(avlevel+i,1)/2 
            av = cat(1,av(1,:,:,:),av,av(end,:,:,:));
        end
        if band_size(avlevel+i-1,2)/4 ~= band_size(avlevel+i,2)/2 
            av = cat(2,av(:,1,:,:),av,av(:,end,:,:));
        end
        if band_size(avlevel+i-1,3)/4 ~= band_size(avlevel+i,3)/2 
            av = cat(3,av(:,:,1,:),av,av(:,:,end,:));
        end

        % Downsample the affine vectors avec by averaging each 2*2*2 group.
        sav = size(av);
        t1 = 1:2:sav(1);
        t2 = 1:2:sav(2);
        t3 = 1:2:sav(3);
        av = 0.125*(av(t1,t2,t3,:) + av(t1,t2+1,t3,:) + av(t1,t2,t3+1,:) + av(t1,t2+1,t3+1,:) + av(t1+1,t2,t3,:) + av(t1+1,t2+1,t3,:) + av(t1+1,t2,t3+1,:) + av(t1+1,t2+1,t3+1,:));
    end
end

% Create x y z indices, going from -1 to +1 across the image.
s1 = size(av,1);
s2 = size(av,2);
s3 = size(av,3);
row_index_distance = 2/band_size(1,1)*(2^lev);
column_index_distance = 2/band_size(1,2)*(2^lev);
frame_index_distance = 2/band_size(1,3)*(2^lev);

xi_temp = ones(s1,1) * [-(s2-1)/2*column_index_distance:column_index_distance:(s2-1)/2*column_index_distance];
yi_temp = [-(s1-1)/2*row_index_distance:row_index_distance:(s1-1)/2*row_index_distance]' * ones(1,s2);

xi = zeros([size(xi_temp),s3]);
yi = zeros([size(yi_temp),s3]);
zi = zeros([size(xi_temp),s3]);
for k = 1:s3,
    xi(:,:,k) = xi_temp;
    yi(:,:,k) = yi_temp;
    zi(1:s1,1:s2,k) = -(s3-1)/2*frame_index_distance + (k-1)*frame_index_distance;
end

%Calc. the elements of shift vector 
shift_y = av(:,:,:,1) + av(:,:,:,4).*yi + av(:,:,:,7).*xi + av(:,:,:,10).*zi;
shift_x = av(:,:,:,2) + av(:,:,:,5).*yi + av(:,:,:,8).*xi + av(:,:,:,11).*zi;
shift_z = av(:,:,:,3) + av(:,:,:,6).*yi + av(:,:,:,9).*xi + av(:,:,:,12).*zi;

return;