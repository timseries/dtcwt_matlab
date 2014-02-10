function dtcwt3dC_xa(level, ext_mode)

% function dtcwt3dC(level)
% Perform the 3-dimensional Dual-Tree Complex Wavelet Transform (DT CWT)
% or its inverse on dataset xa, using quarter-sample orthonormal filters at levels >= 2.
% xa is originally taken as the input to the DT CWT. After the transform,
% xa is the LoLoLo band at the highest level, and ya contains the DTCWT
% coefs at each level.

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
%  global xa h0o h1o g0o g1o h0a h1a g0a g1a
%  load3D_xa('sphere'); for k=1:3; xa=cat(k,xa,0*xa); end % Load a dataset and double its dimensions
%  for k=1:3; dtcwt3dC_xa(k,ext_mode); oct_cplx_xa(k); end   % Do 3 levels of forward transform
%  movi16(xa);  % View the transformed subbands (in complex format)
%  for k=[-3:-1]; oct_cplx_xa(k); dtcwt3dC_xa(k,ext_mode); end   % Do 3 levels of inverse transform
%  sxa=size(xa)/2; xa=xa(1:sxa(1),1:sxa(2),1:sxa(3)); add3D_xa('sphere',-1); max(abs(xa(:)))
% Final answer should be < 1e-11, thus showing perfect reconstruction.

global h0o h1o g0o g1o h0a h1a g0a g1a
global xa   % xa is updated by the LLL_band at every level before next level transform
global ya   % ya is used to store the DTCWT coefs of xa. 
            % ya is in the format of ya{level}{band 1:7}.
            % e.g. ya{3}{1} means the 'HLL' band at level 3.
            % band number: 1 => HLL, 2 => LHL, 3 => HHL, 
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

% sxa = size(xa) / (2^(abs(level)-1));
% if any(rem(sxa,4)),
%    error('Size of xa must be a multiple of 4');
% end
% t1 = 1:sxa(1); t2 = 1:sxa(2); t3 = 1:sxa(3);
% 
% sr = sxa/2;
% s1 = 1:sr(1); s2 = 1:sr(2); s3 = 1:sr(3);
% s1b = s1 + sr(1); s2b = s2 + sr(2); s3b = s3 + sr(3);



fprintf(1,'Level %d: ',level);
tic

if level == 1,
    original_size = size(xa/2); % Original size is the size of the input dataset
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
    
    sxa = size(xa);
    t1 = 1:sxa(1); t2 = 1:sxa(2); t3 = 1:sxa(3);
    sr = sxa/2;
    s1 = 1:sr(1); s2 = 1:sr(2); s3 = 1:sr(3);
    s1b = s1 + sr(1); s2b = s2 + sr(2); s3b = s3 + sr(3);
    % Loop for each slice, incrementing 2nd dimension.
    % (Use the 2nd dim as it is faster to leave the 1st dim as the
    % matrix column index.)
    for f = s2,
        y = reshape(xa(s1,f,s3),sr([1 3])).';
        % Do odd top-level filters on 3rd dim.
        xa(s1,f,s3) = colfilter(y,h0o).';
        xa(s1,f,s3b) = colfilter(y,h1o).';
    end
    
    % Loop for each frame, incrementing 3rd dimension.
    for f = t3,
        % Do odd top-level filters on rows.
        y = xa(s1,s2,f).';
        y2 = [colfilter(y,h0o);  colfilter(y,h1o)].';
        % Do odd top-level filters on columns.
        xa(s1,t2,f) = colfilter(y2,h0o);
        xa(s1b,t2,f) = colfilter(y2,h1o);
    end
    
    for band = 1:7
        % Initialise ya
        ya{level}{band} = zeros(sxa/2); % Size of each of the octal bands must be extended_LLL_band_size/2
    end
    
    % Now record ya{level}{band 1:7} and update LLL_band
    ya{level}{1} = xa(s1,s2b,s3);     % HLL
    ya{level}{2} = xa(s1b,s2,s3);     % LHL
    ya{level}{3} = xa(s1b,s2b,s3);    % HHL
    ya{level}{4} = xa(s1,s2,s3b);     % LLH
    ya{level}{5} = xa(s1,s2b,s3b);    % HLH
    ya{level}{6} = xa(s1b,s2,s3b);    % LHH
    ya{level}{7} = xa(s1b,s2b,s3b);   % HHH
    xa = xa(s1,s2,s3);                % LLL
    
elseif level >= 2,
    % Check if the LoLoLo band is divisable by the value of ext_mode in
    % each direction. If not, we aim to extend the LoLoLo band to make its
    % size divisable of the value of ext_mode.
    LLL_band_size = size(xa);
    if ext_mode == 4
        if any(rem(LLL_band_size(1),4)),	% sr is the size of the LoLoLo band	
           % Extend by 2 rows for each of the octal bands, if no. of rows of LoLoLo is not divisable by 4. 
           xa = cat(1,xa(1,:,:),xa,xa(end,:,:));
        end 
        if any(rem(LLL_band_size(2),4)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 4. 
           xa = cat(2,xa(:,1,:),xa,xa(:,end,:));
        end 
        if any(rem(LLL_band_size(3),4)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 4. 
           xa = cat(3,xa(:,:,1),xa,xa(:,:,end));
        end 
    elseif ext_mode == 8
        if any(rem(LLL_band_size(1),8)),	% sr is the size of the LoLoLo band	
           % Extend by 2 rows for each of the octal bands, if no. of rows of LoLoLo is not divisable by 8. 
           xa = cat(1,xa(1,:,:),xa(1,:,:),xa,xa(end,:,:),xa(end,:,:));
        end 
        if any(rem(LLL_band_size(2),8)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 8. 
           xa = cat(2,xa(:,1,:),xa(:,1,:),xa,xa(:,end,:),xa(:,end,:));
        end 
        if any(rem(LLL_band_size(3),8)),	% sr is the size of the LoLoLo band	
           % Extend by 2 columns for each of the octal bands, if no. of columns of
           % LoLoLo is not divisable by 8. 
           xa = cat(3,xa(:,:,1),xa(:,:,1),xa,xa(:,:,end),xa(:,:,end));
        end
    end
    % Now LLL_band is guaranteed to be divisable by ext_mode
    
    extended_LLL_band_size = size(xa);
    t1 = 1:extended_LLL_band_size(1); t2 = 1:extended_LLL_band_size(2); t3 = 1:extended_LLL_band_size(3);
    s1 = 1:extended_LLL_band_size(1)/2; s2 = 1:extended_LLL_band_size(2)/2; s3 = 1:extended_LLL_band_size(3)/2;
    s1b = s1 + extended_LLL_band_size(1)/2; s2b = s2 + extended_LLL_band_size(2)/2; s3b = s3 + extended_LLL_band_size(3)/2;
    
    for band = 1:7
        % Initialise ya
        ya{level}{band} = zeros(extended_LLL_band_size/2); % Size of each of the octal bands must be extended_LLL_band_size/2
    end
    
    % Loop for each slice, incrementing 2nd dimension.
    for f = t2,
        y = reshape(xa(t1,f,t3),extended_LLL_band_size([1 3])).';
        % Do even Qshift filters on 3rd dim.
        xa(t1,f,s3) = coldfilt(y,h0b,h0a).';
        xa(t1,f,s3b) = coldfilt(y,h1b,h1a).';
    end
    
    % Loop for each frame, incrementing 3rd dimension.
    for f = t3,
        % Do even Qshift filters on rows.
        y = xa(t1,t2,f).';
        y2 = [coldfilt(y,h0b,h0a);  coldfilt(y,h1b,h1a)].';
        % Do even Qshift filters on columns.
        xa(s1,t2,f) = coldfilt(y2,h0b,h0a);
        xa(s1b,t2,f) = coldfilt(y2,h1b,h1a);
    end
    
    % Now record ya{level}{band 1:7} and update LLL_band
    ya{level}{1} = xa(s1,s2b,s3);     % HLL
    ya{level}{2} = xa(s1b,s2,s3);     % LHL
    ya{level}{3} = xa(s1b,s2b,s3);    % HHL
    ya{level}{4} = xa(s1,s2,s3b);     % LLH
    ya{level}{5} = xa(s1,s2b,s3b);    % HLH
    ya{level}{6} = xa(s1b,s2,s3b);    % LHH
    ya{level}{7} = xa(s1b,s2b,s3b);   % HHH
    xa = xa(s1,s2,s3);                % LLL
    
elseif level == -1,    
    sxa = size(xa);
    t1 = 1:2*sxa(1); t2 = 1:2*sxa(2); t3 = 1:2*sxa(3);
    s1 = 1:sxa(1); s2 = 1:sxa(2); s3 = 1:sxa(3);
    s1b = s1 + sxa(1); s2b = s2 + sxa(2); s3b = s3 + sxa(3);
    
    % Now define a temporary variable which is used to combine the highpass
    % bands with the LLL band    
    for k=1:3; xa=cat(k,xa,0*xa); end
    xa(s1,s2b,s3) = ya{abs(level)}{1};
    xa(s1b,s2,s3) = ya{abs(level)}{2};
    xa(s1b,s2b,s3) = ya{abs(level)}{3};
    xa(s1,s2,s3b) = ya{abs(level)}{4};
    xa(s1,s2b,s3b) = ya{abs(level)}{5};
    xa(s1b,s2,s3b) = ya{abs(level)}{6};
    xa(s1b,s2b,s3b) = ya{abs(level)}{7};    
    
    for f = t3,
        % Do odd top-level filters on rows.
        y = colfilter(xa(t1,s2,f).',g0o) + colfilter(xa(t1,s2b,f).',g1o);
        % Do odd top-level filters on columns.
        xa(s1,s2,f) = colfilter(y(:,s1).',g0o) + colfilter(y(:,s1b).',g1o);
    end
    
    for f = s2,
        % Do odd top-level filters on 3rd dim.
        y = reshape(xa(s1,f,t3),sxa(1),2*sxa(3)).';
        xa(s1,f,s3) = (colfilter(y(s3,:),g0o) + colfilter(y(s3b,:),g1o)).';
    end
    
elseif level <= -2,
    extended_LLL_band_size = size(xa);
    t1 = 1:2*extended_LLL_band_size(1); t2 = 1:2*extended_LLL_band_size(2); t3 = 1:2*extended_LLL_band_size(3);
    s1 = 1:extended_LLL_band_size(1); s2 = 1:extended_LLL_band_size(2); s3 = 1:extended_LLL_band_size(3);
    s1b = s1 + extended_LLL_band_size(1); s2b = s2 + extended_LLL_band_size(2); s3b = s3 + extended_LLL_band_size(3);
    
    % Now define a temporary variable which is used to combine the highpass
    % bands with the LLL band
    for k=1:3; xa=cat(k,xa,0*xa); end
    xa(s1,s2b,s3) = ya{abs(level)}{1};
    xa(s1b,s2,s3) = ya{abs(level)}{2};
    xa(s1b,s2b,s3) = ya{abs(level)}{3};
    xa(s1,s2,s3b) = ya{abs(level)}{4};
    xa(s1,s2b,s3b) = ya{abs(level)}{5};
    xa(s1b,s2,s3b) = ya{abs(level)}{6};
    xa(s1b,s2b,s3b) = ya{abs(level)}{7};
    
    for f = t3,
        % Do even Qshift filters on rows.
        y = colifilt(xa(t1,s2,f).',g0b,g0a) + colifilt(xa(t1,s2b,f).',g1b,g1a);
        % Do even Qshift filters on columns.
        xa(t1,t2,f) = colifilt(y(:,s1).',g0b,g0a) + colifilt(y(:,s1b).',g1b,g1a);
    end

    for f = t2,
        y = reshape(xa(t1,f,t3),2*extended_LLL_band_size([1 3])).';
        % Do even Qshift filters on 3rd dim.
        xa(t1,f,t3) = (colifilt(y(s3,:),g0b,g0a) + colifilt(y(s3b,:),g1b,g1a)).';
    end
     
    % Now check if the size of the previous level is exactly twice the size
    % of the current level. If YES, this means we have not done the
    % extension in the previous level. If NO, then we have to remove the
    % appended row / column / frame from the previous level DTCWT coefs.
    size_curr_level = size(ya{abs(level)}{1});   % The band size of the current level
    size_prev_level = size(ya{abs(level)-1}{1}); % The band size of the previous level
    if ext_mode == 4,
        if size_prev_level(1)/size_curr_level(1) ~= 2 
            xa = xa(2:end-1,:,:); % Discard the top and bottom rows
        end
        if size_prev_level(2)/size_curr_level(2) ~= 2
            xa = xa(:,2:end-1,:); % Discard the left and right columns
        end
        if size_prev_level(3)/size_curr_level(3) ~= 2
            xa = xa(:,:,2:end-1); % Discard the left and right frames
        end
    elseif ext_mode == 8,
        if size_prev_level(1)/size_curr_level(1) ~= 2 
            xa = xa(3:end-2,:,:); % Discard the top and bottom rows
        end
        if size_prev_level(2)/size_curr_level(2) ~= 2
            xa = xa(:,3:end-2,:); % Discard the left and right columns
        end
        if size_prev_level(3)/size_curr_level(3) ~= 2
            xa = xa(:,:,3:end-2); % Discard the left and right frames
        end
    end
    
else
   disp('Illegal level');
end

tk = toc;
fprintf(1,'toc = %.2f sec\n',tk);

return;

         
