function [H,g] = ladgen(hh,zthresh)

% function [H,g] = ladgen(hh,zthresh)
% Generate the ladder coefficient matrix H and the gain vector g
% for a ladder implementation of a pair of P-R filters represented
% by the 2 columns of matrix hh.  
%
%  x1 -- g(1,1) + ----------- + --------------- - - --- y1
%      \ g(1,2)     |         |         | 
%       X         H(1,:)    H(2,:)    H(3,:)
%      / g(2,1)     |         |         |
%  x2 -- g(2,2) + - + ----------------- + ----- - - --- y2
%
% x1 is assumed to represent coefs 1,3,5,7... of the input row vector.
% x2 is assumed to represent coefs 2,4,6,8... of the input row vector.
% The 2x2 g matrix specifies how x1 and x2 are mixed into y1 and y2. 
% zthresh is the threshold, below which coefs are regarded as zero.

% Nick Kingsbury, Cambridge University, March 1998.
% 
% To test the results of this routine use:
% [H,g] = ladgen2(beside(h0e,h1e));
% x1=zeros(1,16);x2=x1;x1(8)=1;[y1,y2] = ladfil(x1,x2,H,g,-1);
% y=[y1;y2];beside(y(:),2*g0e)
% x1=zeros(1,16);x2=x1;x2(8)=1;[y1,y2] = ladfil(x1,x2,H,g,-1);
% y=[y1;y2];beside(y(:),2*g1e)
%
% [H,g] = ladgen2([1 2 1 0 0;-1 -2 6 -2 -1]')
% x1=zeros(2,7);x2=x1;x1(1,4)=1;x2(2,4)=1;[y1,y2] = ladfil(x1,x2,H,g);y1,y2

if nargin < 2, zthresh = 1e-6; end

s = findends(hh,zthresh);
Hwidth = 7;
H = zeros(1,Hwidth);
hc0 = (Hwidth+1)/2;
swap = 0; % Defines whether columns of hh have been swapped.
he = 1;  % Defines which end of hh we are reducing (1 = start, 2 = end).

while (any(diff(s) > 0)),  % Continue until h1 and h2 are each length 1.

% If h1 longer than h2, swap columns of hh and increment rows in H.
  if diff(diff(s)) < 0, 
    hh = hh(:,[2 1]);
    H = [zeros(1,Hwidth); H];
    swap = 1 - swap;
    s = findends(hh,zthresh);
  end

% Reduce h2 by subtracting a shifted and scaled version of h1.
  r = hh(s(he,2),2) / hh(s(he,1),1); % Scale factor.
  shift = -diff(s(he,:));  % Shift - should be even if filters are PR.
  if rem(shift,2) ~= 0,
    disp('Warning, shift is odd!');shift % hh
  end
  H(1,hc0-shift/2) = r;
  ta = [s(1,1):s(2,1)];
  tb = ta - shift;
  hh(tb,2) = hh(tb,2) - r * hh(ta,1);
  % hh, H, s, shift, pause
  he = 3 - he;
  s = findends(hh,zthresh);
end

% Add an extra row of zeros to H and swap again if swap is true.
if swap,
  hh = hh(:,[2 1]);
  H = [zeros(1,Hwidth); H];
%  hh, H, pause
end

s = findends(hh,zthresh);
ss = [max(max(s)) min(min(s))];
if diff(ss)~=-1,
    disp('Warning, final hh matrix is not 2x2!');hh
end
g = hh(ss,:).';

% Trim zeros from last row and symmetrical outer columns of H.
sH = size(H);
if all(abs(H(sH(1),:)) < 1e-8), H = H(1:(sH(1)-1),:); end
while all(all(abs(H(:,[1 sH(2)])) < 1e-8)),
  H = H(:,2:(sH(2)-1));
  sH = size(H);
end
return


