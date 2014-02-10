function s = findends(x,zthresh)

% function s = findends(x,zthresh)
% Find the ends of columns of matrix x, defined by the first and last
% non-zero terms in each column.
% s(1,:) is the index of the first non-zero element in each column;
% s(2,:) is the index of the last non-zero element in each column.
% zthresh is the threshold, below which elements are regarded as zero
% and has a default value of 10^(-8).

% Nick Kingsbury, Cambridge University, Feb 98.

sx = size(x);
if nargin < 2, zthresh = 1e-8; end

s = zeros(2,sx(2));

for col=1:sx(2),
  v = x(:,col);
  nzi = find(abs(v) > zthresh);
  s(1,col) = min(nzi);
  s(2,col) = max(nzi);
end

return


