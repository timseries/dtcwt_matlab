% reg3d_V2
% The main 3-D registration programme (Version 2)
% Version 2 of the 3-D registration programme uses DT-CWT to estimate
% the motion affine parameters. Once the motion field is acquired, use it
% to generate the full-resolution motion field and shift the dataset xa in
% spatial domain.

% estshhemm33 performs motion estimation, which generates the affine parameters, based on Nick Kingsbury's algorithm as described in
% "multi-scale displacement estimation and registration for 2-D and 3-D
% datasets", May 2007

% Huizhong Chen & Nick Kingsbury, June 2009

global xa xb ya xa_before_reg ya_before_reg % ya stores the coefs of the highpass bands after xa is DT-CWT transformed
global band_size    % band_size stores the size of the band at each level

% Load and display dataset xa and xb
load3D_xa('NMRset1');
figure(1); settitle('A image'); 
movi16(xa); pause(0.1);
load3D_xb('NMRset5');
figure(2); settitle('B image (reference image)');
movi16(xb); pause(0.1);

figure(3); settitle('A image');
movi16(xa); pause(0.1);

% Set up expected rotation rates for each subband with freq offset.
% To centre the passband of h on the scaling func and wavelet passbands.
jw0 = -pi/2.1;
w = [jw0 jw0*3];

% Specify the parameters and modes for the registration operation
nit = 6;
levelsel = [5 4; 5 4; 5 3; 5 3; 5 3; 5 3];

% mode(1) = 1 means using a global affine model for iteration 1.
% mode(2) = 1 means we assign the constraints of the outer most locality to
% be zero. mode(2) = 2 means we zero the constraints of the subbands
% which are nearly parallel to the row, column and frame direction
% respectively.
mode = [1 2];        % Select the mode of estshhemm33
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

% Clear variables to save memory
clear yb ya ya_before_reg xb
xa = []; xa_before_reg = [];

% Obtain the full-resolution motion field for the registration image
affineshift3_spatial(avec,avlevel);

% Load xa with the dataset to be registered
load3D_xa('NMRset1')
xa_before_reg = xa;              % xa_before_reg will be used to be interpolated to obtain the registered dataset
xa = zeros(size(xa_before_reg)); % initialise xa
spatial_transform;               % Perform 3-D interpolation of xa in spatial domain
% Display the registered image
figure(3); settitle('Registered image');
movi16(xa); pause(0.1);