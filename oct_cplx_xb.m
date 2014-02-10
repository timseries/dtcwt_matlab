function oct_cplx_xb(level)
% function oct_cplx_xa(level)
% Convert between 3-D octets and four subbands of complex coefs at a given
% level of the 3-D DT CWT.
% When in 4-subband complex form, the 4 subbands are spatially separated
% within each frame of xb, and adjacent frames (along 3rd dim) are the real
% and imag parts of each complex sample.
%
% Forward transforms (convert octets to 4 complex subbands):
%  level > 0
%
% Inverse transforms (convert 4 complex subands to octets):
%  level < 0
%
% Functions o2c() and c2o() are used for converting each subband in
% the appropriate direction.

% Nick Kingsbury, Cambridge University, July 1999.
% Modified by Huizhong Chen, Oct 2008

if level > 0,  % Convert from octet to complex.
    o2c(level); % o2c converts the 7 highpass bands from octet to complex    
else  % Convert from complex to octet.
    c2o(level);
end
return

function o2c(level)
% function o2c(level)
% Convert the 7 bands of octets numbers to complex for current level.
global yb

sya = size(yb{level}{1});
t3 = 1:2:sya(3);
s1 = 1:sya(1); s2 = 1:sya(2);
for band = 1:7
    for frame = t3
        f2 = frame + [0 1];
        y1 = squeeze(yb{level}{band}(s1,s2,f2(1)));
        y2 = squeeze(yb{level}{band}(s1,s2,f2(2)));

        sy = size(y1);
        t1 = 1:2:sy(1); t2 = 1:2:sy(2);

        % Arrange voxels from the corners of the octets into
        % 8 separate subimages.
        %  a----b
        %  |\   |\
        %  | e----f
        %  c-|--d |
        %   \|   \|
        %    g----h
        a = y1(t1,t2);
        b = y1(t1,t2+1);
        c = y1(t1+1,t2);
        d = y1(t1+1,t2+1);
        e = y2(t1,t2);
        f = y2(t1,t2+1);
        g = y2(t1+1,t2);
        h = y2(t1+1,t2+1);

        % Form the 4 subbands of real parts.
        yb{level}{band}(s1,s2,f2(1)) = [a-d-f-g  a+d+f-g; a+d-f+g  a-d+f+g]*0.5;
        % Form the 4 subbands of imag parts.
        yb{level}{band}(s1,s2,f2(2)) = [b+c+e-h  -b+c+e+h; b-c+e+h  -b-c+e-h]*0.5;
    end
end
return

function c2o(level)
% function c2o(level)
% Convert the 7 bands of complex numbers to octets for current level.
global yb
sya = size(yb{abs(level)}{1});
t3 = 1:2:sya(3);
s1 = 1:sya(1); s2 = 1:sya(2);
for band = 1:7
    for frame = t3
        f2 = frame + [0 1];
        y1 = squeeze(yb{abs(level)}{band}(s1,s2,f2(1)));  % Real
        y2 = squeeze(yb{abs(level)}{band}(s1,s2,f2(2)));  % Imag

        sy = size(y1)/2;
        t1 = 1:sy(1); t2 = 1:sy(2);

        % Arrange voxels from the real and imag parts of the 4 subbands
        % into 8 separate subimages .
        %  A----B   Real layer
        %  |\   |\
        %  | E----F   Imag layer
        %  C-|--D |
        %   \|   \|
        %    G----H
        A = y1(t1,t2);
        B = y1(t1,t2+sy(2));
        C = y1(t1+sy(1),t2);
        D = y1(t1+sy(1),t2+sy(2));
        E = y2(t1,t2);
        F = y2(t1,t2+sy(2));
        G = y2(t1+sy(1),t2);
        H = y2(t1+sy(1),t2+sy(2));

        t1 = s1(1:2:length(s1));
        t2 = s2(1:2:length(s2));

        % Recover each of the 8 corners of the octets.
        yb{abs(level)}{band}(t1,t2,f2(1))     =  (A+B+C+D)*0.5; % a 
        yb{abs(level)}{band}(t1,t2+1,f2(1))   =  (E-F+G-H)*0.5; % b 
        yb{abs(level)}{band}(t1+1,t2,f2(1))   =  (E+F-G-H)*0.5; % c
        yb{abs(level)}{band}(t1+1,t2+1,f2(1)) = (-A+B+C-D)*0.5; % d

        yb{abs(level)}{band}(t1,t2,f2(2))     =  (E+F+G+H)*0.5; % e 
        yb{abs(level)}{band}(t1,t2+1,f2(2))   = (-A+B-C+D)*0.5; % f 
        yb{abs(level)}{band}(t1+1,t2,f2(2))   = (-A-B+C+D)*0.5; % g
        yb{abs(level)}{band}(t1+1,t2+1,f2(2)) = (-E+F+G-H)*0.5; % h
    end
end
return
