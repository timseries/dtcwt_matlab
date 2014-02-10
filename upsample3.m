function y = upsample3(x,factor,h1,h2,h3)

% function y = upsample3(x,factor,h1,h2,h3)
%
% Up-sample the matrix x by factor(1:3) in each direction.
% h1, h2 and h3, if given, are the filtering functions for columns and rows and frames.
% factor defaults to 2*2*2 up sampling, if not given.
% h1 h2 and h3 default to bilinear interpolation. 
%
% Huizhong Chen, Cambridge University, July 2008.

% Default to 2:1 upsampling.
if nargin < 2, factor = [2 2 2]; end
% Default to equal upsampling in 3 directions if factor is a scalar.
if length(factor) < 2, factor = factor*[1 1 1]; end

if all(factor == [1 1 1])
   y = x;
   return;
end

if nargin < 3,  % Define h1 h2 and h3 as rectangular pulses for bilinear interpolation.
   h1 = ones(factor(1),1) / factor(1);
   h2 = ones(factor(2),1) / factor(2);
   h3 = ones(factor(3),1) / factor(3);
   % Adjust to make lengths of h1 h2 and h3 odd to get correct alignment.
   if rem(factor(1),2) == 0, h1 = conv(h1,[1;1]/2); end
   if rem(factor(2),2) == 0, h2 = conv(h2,[1;1]/2); end
   if rem(factor(3),2) == 0, h3 = conv(h3,[1;1]/2); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sx = size(x);
s1 = sx(1); s2 = sx(2); s3=sx(3);
y = zeros(s1*factor(1),s2*factor(2),s3*factor(3));
for k=1:s3
    for i = 1:factor(3)
        y(:,:,(k-1)*factor(3)+i) = kron(x(:,:,k),ones(factor(1:2)));
    end
end %Now y has been upsampled in 3 dims by a factor of 'factor'

%Filter the y matrix
sy = size(y);
s1 = 1:sy(1); s2 = 1:sy(2); s3 = 1:sy(3);
for k = s3,
    z = y(:,:,k).';
    y(s1,s2,k) = colfilter(z,h1).'; %filtering on rows
    y(s1,s2,k) = colfilter(y(s1,s2,k),h2); %filtering on columns
end

for k = s2,
    z = reshape(y(s1,k,s3),sy([1 3])).';
    z = colfilter(z,h3); %filtering on columns

    y(s1,k,s3) = z(s3,:)';
end
return


