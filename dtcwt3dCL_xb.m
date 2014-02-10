function dtcwt3dCL_xb(level,ext_mode)

% function dtcwt3dCL_xb(level,ext_mode)
% Perform the 3-dimensional Dual-Tree Complex Wavelet Transform (DT CWT)
% or its inverse on dataset xb, using quarter-sample orthonormal filters at levels >= 2.
% xb is originally taken as the input to the DT CWT. After the transform,
% xb is the LoLoLo band at the highest level, and yb contains the DTCWT
% coefs at each level.
% NOTE: For efficiency this version ignores all high bands at level 1.

%
% Forward transforms:
%  Top level: level = 1
%  Levels >= 2: level >= 2
%
% Inverse transforms:
%  Top level: level = -1
%  Levels >= 2: level <= -2
% 
% There are two values for ext_mode, either 4 or 8.
% If ext_mode = 4, check whether 1st level is divisible by 2 (if not we
% give an error msg). Also check whether from 2nd level onwards, the coefs can
% be divided by 4. If any dimension size is not a multiple of 4, append
% extra coefs by repeating the edges.
% If ext_mode = 8, check whether 1st level is divisible by 4 (if not we
% give an error msg). Also check whether from 2nd level onwards, the coefs can
% be divided by 8. If any dimension size is not a multiple of 8, append
% extra coefs by repeating the edges twice.
%
% Nick Kingsbury, Cambridge University, July 1999.
% This version modified to use column filtering, Nov 2003.
% ext_mode added by Huizhong Chen, Jan 2009

% To test use:
%  global xb h0o h1o g0o g1o h0a h1a g0a g1a
%  load3d_xb('sphere'); % Load a dataset.
%  dtcwt3dCL_xb(1, ext_mode); for k=2:3; dtcwt3dCL_xb(k,ext_mode); oct_cplx_xb(k); end   % Do 3 levels of forward transform
%  for k= -3:-2; oct_cplx_xb(k); dtcwt3dCL_xb(k,ext_mode); end; dtcwt3dCL_xb(-1,ext_mode); % Do 3 levels of inverse transform
%  add3D_xb('sphere',-1); mean(abs(xb(:)))
% Final answer should be < 2, thus showing near-perfect reconstruction, apart from loss of
% level 1 high band components.

global h0o h1o g0o g1o h0a h1a g0a g1a
global xb   % xb is updated by the LLL_band at every level before next level transform
global yb   % yb is used to store the DTCWT coefs of xb. 
            % yb is in the format of yb{level}{band 1:7}, level starts
            % from 2 because the 1st level has only LLL_band.
            % e.g. yb{3}{1} means the 'HLL' band at level 3.
            % 1 => HLL, 2 => LHL, 3 => HHL, 
            % 4 => LLH, 5 => HLH, 6 => LHH, 7 => HHH
global original_size  % original_size is the original input image size

if nargin < 2, ext_mode = 4; end % default ext_mode is 4

% Check whether ext_mode has been given a correct value
if ext_mode~=4 && ext_mode~=8
    error('In dtcwt3dCL, ext_mode must be either 4 or 8');
end

% Define the top level filters.
if isempty(h0o),
   disp('Generating level 1 filters.')
   [h0o,h1o,g0o,g1o] = wavegen('near-sym',[3.5 3/16 0]);
%   [h0o,h1o,g0o,g1o] = wavegen('Antonini');
%   [h0o,h1o,g0o,g1o] = wavegen('LeGall');
end

% Define the main orthogonal Qshift filters.
if isempty(h0a),
   disp('Generating level 2+ filters.')
   [h0a,h1a,g0a,g1a] = wavegen('Qshift',[14 0.36]);
%   [h0a,h1a,g0a,g1a] = wavegen('Qshift',8);
   % These filters have unity gain at zero or fs/2.
   
   % Normalise analysis filter gains to preserve energy in xfm and give PR.
   hscale = 1 / sqrt(sum(h0a .* h0a));
   h0a = h0a * hscale;
   h1a = h1a * hscale;
   g0a = g0a * hscale;
   g1a = g1a * hscale;
end

h0b = rev(h0a);
h1b = rev(h1a);
g0b = rev(g0a);
g1b = rev(g1a);

fprintf(1,'Level %d: ',level);
tic

if level == 1,
    original_size = size(xb);
    if ext_mode == 4
        if any(rem(original_size,2))
            warning('Image size should have even samples in each direction');  %#ok<WNTAG>
            return;
        end
    elseif ext_mode ==8
        if any(rem(original_size,4))
            warning('Image size should be divisible by 4 in each direction');  %#ok<WNTAG>
            return;
        end
    end

    sxb = size(xb);
    sxb = 2*sxb; % Compensate for lack of high bands at level 1.

    sr = sxb/2;
    s1 = 1:sr(1); s2 = 1:sr(2); s3 = 1:sr(3);

    % Loop for each slice, incrementing 2nd dimension.
    % (Use the 2nd dim as it is faster to leave the 1st dim as the
    % matrix column index.)
    for f = s2,
        y = reshape(xb(s1,f,s3),sr([1 3])).';
        % Do odd top-level filters on 3rd dim.
        xb(s1,f,s3) = colfilter(y,h0o).';
    end
    
    % Loop for each frame, incrementing 3rd dimension.
    for f = s3,
        % Do odd top-level filters on rows.
        y2 = colfilter(xb(s1,s2,f).',h0o).';
        % Do odd top-level filters on columns.
        xb(s1,s2,f) = colfilter(y2,h0o);
    end
    
elseif level >= 2,
    % Check if the LoLoLo band is divisable by 4 in each direction. If
    % not, we aim to extend the LoLoLo band to make its size divisable by
    % 4. 
    LLL_band_size = size(xb);
    if ext_mode == 4
        if any(rem(LLL_band_size(1),4)),	% sr is the size of the LoLoLo band	
           % Extend by 2 rows for each of the octal bands, if no. of rows of LoLoLo is not divisable by 4. 
           xb = cat(1,xb(1,:,:),xb,xb(end,:,:));
        end 
        if any(rem(LLL_band_size(2),4)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 4. 
           xb = cat(2,xb(:,1,:),xb,xb(:,end,:));
        end 
        if any(rem(LLL_band_size(3),4)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 4. 
           xb = cat(3,xb(:,:,1),xb,xb(:,:,end));
        end 
    elseif ext_mode == 8
        if any(rem(LLL_band_size(1),8)),	% sr is the size of the LoLoLo band	
           % Extend by 2 rows for each of the octal bands, if no. of rows of LoLoLo is not divisable by 8. 
           xb = cat(1,xb(1,:,:),xb(1,:,:),xb,xb(end,:,:),xb(end,:,:));
        end 
        if any(rem(LLL_band_size(2),8)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 8. 
           xb = cat(2,xb(:,1,:),xb(:,1,:),xb,xb(:,end,:),xb(:,end,:));
        end 
        if any(rem(LLL_band_size(3),8)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 8. 
           xb = cat(3,xb(:,:,1),xb(:,:,1),xb,xb(:,:,end),xb(:,:,end));
        end
    end
    % Now LLL_band is guaranteed to be divisable by ext_mode
    
    extended_LLL_band_size = size(xb);
    t1 = 1:extended_LLL_band_size(1); t2 = 1:extended_LLL_band_size(2); t3 = 1:extended_LLL_band_size(3);
    s1 = 1:extended_LLL_band_size(1)/2; s2 = 1:extended_LLL_band_size(2)/2; s3 = 1:extended_LLL_band_size(3)/2;
    s1b = s1 + extended_LLL_band_size(1)/2; s2b = s2 + extended_LLL_band_size(2)/2; s3b = s3 + extended_LLL_band_size(3)/2;
    
    for band = 1:7
        % Initialise yb
        yb{level}{band} = zeros(extended_LLL_band_size/2); % Size of each of the octal bands must be extended_LLL_band_size/2
    end
    % Loop for each slice, incrementing 2nd dimension.
    for f = t2,
        y = reshape(xb(t1,f,t3),extended_LLL_band_size([1 3])).';
        % Do even Qshift filters on 3rd dim.  
        xb(t1,f,s3) = coldfilt(y,h0b,h0a).';
        xb(t1,f,s3b) = coldfilt(y,h1b,h1a).';
    end
    
    % Loop for each frame, incrementing 3rd dimension.
    for f = t3,
        % Do even Qshift filters on rows.
        y = xb(t1,t2,f).';
        y2 = [coldfilt(y,h0b,h0a);  coldfilt(y,h1b,h1a)].';
        % Do even Qshift filters on columns.
        xb(s1,t2,f) = coldfilt(y2,h0b,h0a);
        xb(s1b,t2,f) = coldfilt(y2,h1b,h1a);
    end
    
    % Now record yb{level}{band 1:7} and update LLL_band
    yb{level}{1} = xb(s1,s2b,s3);     % HLL
    yb{level}{2} = xb(s1b,s2,s3);     % LHL
    yb{level}{3} = xb(s1b,s2b,s3);    % HHL
    yb{level}{4} = xb(s1,s2,s3b);     % LLH
    yb{level}{5} = xb(s1,s2b,s3b);    % HLH
    yb{level}{6} = xb(s1b,s2,s3b);    % LHH
    yb{level}{7} = xb(s1b,s2b,s3b);   % HHH
    xb = xb(s1,s2,s3);                % LLL
    
elseif level == -1,
    extended_LLL_band_size = size(xb);
    s1 = 1:extended_LLL_band_size(1); s2 = 1:extended_LLL_band_size(2); s3 = 1:extended_LLL_band_size(3);
       
    for f = s3,
        % Do odd top-level filters on rows.
        y = colfilter(xb(s1,s2,f).',g0o);
        % Do odd top-level filters on columns.
        xb(s1,s2,f) = colfilter(y(:,s1).',g0o);
    end
    
    for f = s2,
        % Do odd top-level filters on 3rd dim.
        y = reshape(xb(s1,f,s3),extended_LLL_band_size(1),extended_LLL_band_size(3)).';
        xb(s1,f,s3) = colfilter(y(s3,:),g0o).';
    end
    
elseif level <= -2,
    extended_LLL_band_size = size(xb);
    t1 = 1:2*extended_LLL_band_size(1); t2 = 1:2*extended_LLL_band_size(2); t3 = 1:2*extended_LLL_band_size(3);
    s1 = 1:extended_LLL_band_size(1); s2 = 1:extended_LLL_band_size(2); s3 = 1:extended_LLL_band_size(3);
    s1b = s1 + extended_LLL_band_size(1); s2b = s2 + extended_LLL_band_size(2); s3b = s3 + extended_LLL_band_size(3);
    
    % Now define a temporary variable which is used to combine the highpass
    % bands with the LLL band
    for k=1:3; xb=cat(k,xb,0*xb); end
    xb(s1,s2b,s3) = yb{abs(level)}{1};
    xb(s1b,s2,s3) = yb{abs(level)}{2};
    xb(s1b,s2b,s3) = yb{abs(level)}{3};
    xb(s1,s2,s3b) = yb{abs(level)}{4};
    xb(s1,s2b,s3b) = yb{abs(level)}{5};
    xb(s1b,s2,s3b) = yb{abs(level)}{6};
    xb(s1b,s2b,s3b) = yb{abs(level)}{7};

    for f = t3,
        % Do even Qshift filters on rows.
        y = colifilt(xb(t1,s2,f).',g0b,g0a) + colifilt(xb(t1,s2b,f).',g1b,g1a);
        % Do even Qshift filters on columns.
        xb(t1,t2,f) = colifilt(y(:,s1).',g0b,g0a) + colifilt(y(:,s1b).',g1b,g1a);
    end

    for f = t2,
        y = reshape(xb(t1,f,t3),2*extended_LLL_band_size([1 3])).';
        % Do even Qshift filters on 3rd dim.
        xb(t1,f,t3) = (colifilt(y(s3,:),g0b,g0a) + colifilt(y(s3b,:),g1b,g1a)).';
    end
    
    % Now check if the size of the previous level is exbctly twice the size
    % of the current level. If YES, this means we have not done the
    % extension in the previous level. If NO, then we have to remove the
    % appended row / column / frame from the previous level DTCWT coefs.
    size_curr_level = size(yb{abs(level)}{1});
    if(level<=-3)
        size_prev_level = size(yb{abs(level)-1}{1});
        if ext_mode == 4,
            if size_prev_level(1)/size_curr_level(1) ~= 2 
                xb = xb(2:end-1,:,:); % Disgard the top and bottom rows
            end
            if size_prev_level(2)/size_curr_level(2) ~= 2
                xb = xb(:,2:end-1,:); % Disgard the left and right columns
            end
            if size_prev_level(3)/size_curr_level(3) ~= 2
                xb = xb(:,:,2:end-1); % Disgard the left and right columns
            end
        elseif ext_mode == 8,
            if size_prev_level(1)/size_curr_level(1) ~= 2 
                xb = xb(3:end-2,:,:); % Disgard the top and bottom rows
            end
            if size_prev_level(2)/size_curr_level(2) ~= 2
                xb = xb(:,3:end-2,:); % Disgard the left and right columns
            end
            if size_prev_level(3)/size_curr_level(3) ~= 2
                xb = xb(:,:,3:end-2); % Disgard the left and right columns
            end
        end
    else %when level=-2
        if ext_mode == 4,
            if original_size(1) / size_curr_level(1) ~= 2
                xb = xb(2:end-1,:,:); % Disgard the top and bottom rows
            end
            if original_size(2)/size_curr_level(2) ~= 2
                xb = xb(:,2:end-1,:); % Disgard the left and right columns
            end
            if original_size(3)/size_curr_level(3) ~= 2
                xb = xb(:,:,2:end-1); % Disgard the top and bottom frames
            end
        elseif ext_mode == 8,
            if original_size(1) / size_curr_level(1) ~= 2
                xb = xb(3:end-2,:,:); % Disgard the top and bottom rows
            end
            if original_size(2)/size_curr_level(2) ~= 2
                xb = xb(:,3:end-2,:); % Disgard the left and right columns
            end
            if original_size(3)/size_curr_level(3) ~= 2
                xb = xb(:,:,3:end-2); % Disgard the top and bottom frames
            end
        end
    end
    

else
   disp('Illegal level');
end

tk = toc;
fprintf(1,'toc = %.2f sec\n',tk);

return;