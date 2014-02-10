% reg3d_V1
% The main 3-D registration programme (Version 1)
% Version 1 of the 3-D registration programme uses DT-CWT to estimate
% the motion affine parameters. Once the motion field is acquired, use it
% to shift the DT-CWT subband coefficients and then inverse DT-CWT to register
% dataset xa to xb.

% estshhemm33 performs motion estimation, which generates the affine parameters, based on Nick Kingsbury's algorithm as described in
% "multi-scale displacement estimation and registration for 2-D and 3-D
% datasets", May 2007

% Huizhong Chen & Nick Kingsbury, June 2009

global xa xb ya ya_before_reg % ya stores the coefs of the highpass bands after xa is DT-CWT transformed
global band_size    % band_size stores the size of the band at each level

% Load and display dataset xa and xb
load3D_xa('sphere_shifted');
figure(1); settitle('A image'); 
movi16(xa); pause(0.1);
load3D_xb('sphere');
figure(2); settitle('B image (reference image)');
movi16(xb); pause(0.1);

figure(3); settitle('A image');
movi16(xa); pause(0.1);

% Set up expected rotation rates for each subband with freq offset.
% To centre the passband of h on the scaling func and wavelet passbands.
jw0 = -pi/2.1;
w = [jw0 jw0*3];

% Specify the parameters and modes for the registration operation
%nit = 6;
%levelsel = [5 4; 5 4; 5 3; 4 3; 4 3; 4 3];
nit = 2;
levelsel = [4 3; 3 2];

% mode(1) = 1 means using a global affine model for iteration 1.
% mode(2) = 1 means we assign the constraints of the outer most locality to
% be zero. mode(2) = 2 means we zero the constraints of the subbands
% which are nearly parallel to the row, column and frame direction
% respectively.
mode = [1 2];        % Select the mode of esthemm33
qscale = 1;          % The scale at which the motion field is plotted
avlevel = 3;         % The affine parameters are to be generated at the scale of the level 'avlevel'
sizeqfilt = [0 4 2]; % User can change size of qfilt as desired

if (size(xa)~=size(xb))
    error('Input images must have same size')
end

max_level = max(max(levelsel));

if any(rem(size(xa),2^max_level))
    ext_mode = 8;    % 8 means the LoLoLo band at each level is to be appended to be a multiple of 8
else
    ext_mode = 4;    % 4 means the LoLoLo band at each level is to be appended to be a multiple of 4
end

band_size = calcbandsize(max_level, ext_mode); % Calculate the octal band size at each level

% Now perfrom DTCWT onto xa and xb
dtcwt3dCL_xa(1,ext_mode);
for k=2:max_level
    dtcwt3dCL_xa(k,ext_mode);
    oct_cplx_xa(k);
end   % Do 'max_lev' levels of forward transform on xa

ya_before_reg = ya;

dtcwt3dCL_xb(1,ext_mode);
for k=2:max_level
    dtcwt3dCL_xb(k,ext_mode);
    oct_cplx_xb(k);
end   % Do 'max_lev' levels of forward transform on xb

% Acquire avec which describes the motion from xa to xb
avec = estshhemm33(w, nit, levelsel, avlevel, mode, qscale, sizeqfilt);

% Reload xa with the DTCWT coefs before it is registered
ya = ya_before_reg;

% Shift xa subbands at each level according to avec
for lev = max_level:-1:2
 srefh = band_size(avlevel,:); %srefh contains the size of the sub-band at level = avlevel
 s2 = srefh(1)/2; 
 s1 = srefh(2)/2;
 s3 = srefh(3);
    [shift_x, shift_y, shift_z] = affineshift3(avec,lev,avlevel);
  
      method = 'spline';       
      yi_distance = 2/band_size(1,1)*(2^lev);
      xi_distance = 2/band_size(1,2)*(2^lev);
      zi_distance = 2/band_size(1,3)*(2^lev);
      % Use shift vectors (adjusted for pel units instead of half-image units) to shift subband.
      shift_xa_cwt_bands3('hi', lev,shift_x/xi_distance,shift_y/yi_distance,shift_z/zi_distance,method); % This line shifts the 28 subbands at the specified level
  
       % Shift xal if at 2nd coarest resolution
       if lev == max_level-1,
          shift_xa_cwt_bands3('lo', max_level,shift_x/xi_distance,shift_y/yi_distance,shift_z/zi_distance,method);
       end
end


for k = -max_level:-2,
    oct_cplx_xa(k); 
    dtcwt3dCL_xa(k,ext_mode);
end; 
dtcwt3dCL_xa(-1,ext_mode); % Do max_level levels of inverse transform

% Display the registered image
figure(3); settitle('Registered image');
movi16(xa); pause(0.1);

