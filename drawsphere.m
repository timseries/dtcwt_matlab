function drawsphere(r,w,dx,N)
% function drawsphere(r,w,dx,N)
% Generate a volume of size N(1)*N(2)*N(3) pels, containing a sphere
% radius r*N and centred at dx(1:3) relative
% to the centre of the image.  The edge of the sphere is a cosine shaped
% edge of width w (from 10 to 90% points).
%
% eg: 
% drawsphere(0.4,0.01,0,128);
% generates a sphere of radius 0.4*128 pels and with an edge width of
% 0.01*128 pels.
%
% Nick Kingsbury, Cambridge University, Nov 2003.

if length(N) == 1, N = [1 1 1]*N; end
if length(dx) == 1, dx = [1 1 1]*dx; end

global x

x = zeros(N);

% define pic in each frame and assemble to form a 3-d sphere
for k = 1:N(3),
   ix = ones(N(1),1) * (([1:N(2)] - (N(2)+1)/2 - dx(1))/(r*N(2)));
   iy = (([1:N(1)]' - (N(1)+1)/2 - dx(2))/(r*N(1))) * ones(1,N(2));
   iz = (k - (N(3)+1)/2 - dx(3))/(r*N(3));
   x(:,:,k) = sin(min(max((exp(-0.5 * (ix .* ix + iy .* iy + iz * iz)) - exp(-0.5))*(r*3/w), -pi/2), pi/2));
   x(:,:,k) = x(:,:,k) + 1;
end

return

