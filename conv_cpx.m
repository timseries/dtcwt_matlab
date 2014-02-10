function y = conv_cpx( x )
% y = conv_cpx( x )
% This function converts x (which has real & imag parts of coefs stored in ajacent layer)
% to truely complex array. Therefore the dimension in z direction is reduced by a factor of 2.

sx = size(x);

if rem(sx(3),2)
    error('z direction must have even numbered size');
end

% Initialise y array
y = zeros(sx(1),sx(2),sx(3)/2);

% Now calc. y
for k = 1:sx(3)/2
    y(:,:,k) = x(:,:,2*(k-1)+1) + sqrt(-1)*x(:,:,2*k);
end

return