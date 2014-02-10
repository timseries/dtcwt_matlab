function  band_size = calcbandsize(max_transform_lev, ext_mode)
% function band_size = calcsubbandsize( max_transform_lev, ext_mode )
% Given the original image xa, calculate the size of the band at each
% level. N.B. band is the size of one of the octal bands at the level.
% band_size is a vector with length 'max_transform_lev', and each element
% stores the band size at the level specified by the vector index.
% ext_mode is used to decide how to do extend the size of the signal in
% order to make it a multiple of ext_mode at the LoLoLo band of each level

global xa

sxa = size(xa);

if any(rem(sxa,ext_mode/2))
    error('The input image size must be a multiple of %d in each dimension', ext_mode/2);
end

% initialise
band_size = zeros(max_transform_lev, 3);


band_size(1,:) = sxa; % For level 1 transform, the band_size is equal to the original image size

for lev = 2:max_transform_lev
    sxa = (sxa + rem(sxa,ext_mode))/2;
    band_size(lev,:) = sxa;
end