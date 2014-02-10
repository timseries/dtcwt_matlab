function y = cpxinterp2b(x,locs,w,interpmethod)

% function y = cpxinterp2b(x,locs,w,interpmethod)
%
% Complex bandpass interpolation routine for the matrix x.
% Uses bandpass interpolation centred on w(1:2) rad/sample centre freq.
% locs(:,1:2) specifies the locations of the interpolated points
% in the matrix x. (eg locs(1,:) specifies the row and column locations of the first sample.)
% w specifies the approximate expected rotation in radians between consecutive samples 
% within the matrix x.  w(1) is the rate of rotation down the columns, and w(2) 
% across the rows of x, so that phase differences can be unwrapped correctly.
% interpmethod is a string which specifies the rule used by interp2:
%   'nearest','linear','cubic','spline'  are valid methods in
% increasing order of computational complexity (defaults to 'cubic').
% 'linear' is about 4 times as fast as 'cubic'.
%
% Version 2b is more similar to the function interp2, in that it lists arbitrary locations
% for the interpolated points rather than a regular upsampling pattern.
%
% For the DTCWT, the w matrix for the 6 subbands should be:
%  w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally pi/2, but reduced a bit due to asymmetry of subband freq responses.
%
% Nick Kingsbury, Cambridge University, August 2005

if nargin < 4, interpmethod = 'cubic'; end

nr = size(x,1);
nc = size(x,2);

sl = size(locs);
if sl(2) ~=2, error('cpxinter2b - locs must be 2 columns wide'); end
jw = sqrt(-1) * w;

% Set up matrix padding vectors 1 row/column wide.
z1 = zeros(1,nc+2);
z2 = zeros(nr,1);

% Create linear ramp matrices for phase wrapping.
thx = ones(nr,1) * [1:nc];
thy = [1:nr]' * ones(1,nc);

% Get the image, unwrap the phases and pad it with zeros.
ye = [z1; z2  x.*exp(-thx*jw(2) - thy*jw(1))  z2; z1];

% Interpolate ye to the new points, specified in locs(:,1:2).
yi = interp2(0:(nc+1),0:(nr+1),ye,locs(:,2),locs(:,1),interpmethod);  % Add one to xs,ys to allow for extension of ye.
% Rewrap the phases.
y = yi .* exp(locs(:,2)*jw(2) + locs(:,1)*jw(1));

return;
