function y = abs3D(x,mode)

% function abs3D(x,mode)
% Find the absolute magnitude of x, where real and imag parts are stored
% either as alternate pels along each row of each slice (mode = 2), or in
% 3D mode (mode = 3) as alternate slices.
%
% Nick Kingsbury, Oct 99.


sx = size(x);

if nargin < 2, mode = length(sx); end

t = 1:2:sx(2);
if mode == 2, sx(2) = sx(2)/2; 
else sx(3) = sx(3)/2;
end

y = zeros(sx);

if length(sx) == 3,
	for f = 1:sx(3),
   	if mode == 2,
      	y(:,:,f) = abs(x(:,t,f) + j*x(:,t+1,f));
	   else
   	   y(:,:,f) = abs(x(:,:,2*f-1) + j*x(:,:,2*f));
      end
   end
else
  	y = abs(x(:,t) + j*x(:,t+1));
end

return;
