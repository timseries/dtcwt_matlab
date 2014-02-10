function y = conv_coefs( x )
% y = conv_coefs( x )
% This function converts the true complex element arrays x into real and
% imag. parts, stored in adjacent layers in z direction.
% Therefore the size of z will be doubled.

sx = size(x);

% Initialise y array
y = zeros(sx(1),sx(2),sx(3)*2);

% Now calc. y
for k = 1:sx(3)
    y(:,:,2*k-1) = real(x(:,:,k));
    y(:,:,2*k) = imag(x(:,:,k));
end

return