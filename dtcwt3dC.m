function dtcwt3dC(level)

% function dtcwt3dC(level)
% Perform the 3-dimensional Dual-Tree Complex Wavelet Transform (DT CWT)
% or its inverse, using quarter-sample orthonormal filters at levels >= 2.
% The array of input voxels and of output results is global variable x.
% x must be twice the size (in each dimension) of the real spatial domain data.
%
% Forward transforms:
%  Top level: level = 1
%  Levels >= 2: level >= 2
%
% Inverse transforms:
%  Top level: level = -1
%  Levels >= 2: level <= -2
%
% Nick Kingsbury, Cambridge University, July 1999.
% This version modified to use column filtering, Nov 2003.

% To test use:
%  global x h0o h1o g0o g1o h0a h1a g0a g1a
%  load3d('sphere'); for k=1:3; x=cat(k,x,0*x); end % Load a dataset and double its dimensions
%  for k=1:3; dtcwt3dC(k); oct_cplx(k,0); end   % Do 3 levels of forward transform
%  movi16(x);  % View the transformed subbands (in complex format)
%  for k=[-3:-1]; oct_cplx(k,0); dtcwt3dC(k); end   % Do 3 levels of inverse transform
%  sx=size(x)/2; x=x(1:sx(1),1:sx(2),1:sx(3)); add3d('sphere',-1); max(abs(x(:)))
% Final answer should be < 1e-11, thus showing perfect reconstruction.

global x h0o h1o g0o g1o h0a h1a g0a g1a

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

sx = size(x) / (2^(abs(level)-1));
if any(rem(sx,4)),
   error('Size of x must be a multiple of 4');
end
t1 = 1:sx(1); t2 = 1:sx(2); t3 = 1:sx(3);

sr = sx/2;
s1 = 1:sr(1); s2 = 1:sr(2); s3 = 1:sr(3);
s1b = s1 + sr(1); s2b = s2 + sr(2); s3b = s3 + sr(3);

fprintf(1,'Level %d: ',level);
tic

if level == 1,
        
    % Loop for each slice, incrementing 2nd dimension.
    % (Use the 2nd dim as it is faster to leave the 1st dim as the
    % matrix column index.)
    for f = s2,
        y = reshape(x(s1,f,s3),sr([1 3])).';
        % Do odd top-level filters on 3rd dim.
        x(s1,f,s3) = colfilter(y,h0o).';
        x(s1,f,s3b) = colfilter(y,h1o).';
    end
    
    % Loop for each frame, incrementing 3rd dimension.
    for f = t3,
        % Do odd top-level filters on rows.
        y = x(s1,s2,f).';
        y2 = [colfilter(y,h0o);  colfilter(y,h1o)].';
        % Do odd top-level filters on columns.
        x(s1,t2,f) = colfilter(y2,h0o);
        x(s1b,t2,f) = colfilter(y2,h1o);
    end
    
elseif level >= 2,
    
    % Loop for each slice, incrementing 2nd dimension.
    for f = t2,
        y = reshape(x(t1,f,t3),sx([1 3])).';
        % Do even Qshift filters on 3rd dim.
        x(t1,f,s3) = coldfilt(y,h0b,h0a).';
        x(t1,f,s3b) = coldfilt(y,h1b,h1a).';
    end
    
    % Loop for each frame, incrementing 3rd dimension.
    for f = t3,
        % Do even Qshift filters on rows.
        y = x(t1,t2,f).';
        y2 = [coldfilt(y,h0b,h0a);  coldfilt(y,h1b,h1a)].';
        % Do even Qshift filters on columns.
        x(s1,t2,f) = coldfilt(y2,h0b,h0a);
        x(s1b,t2,f) = coldfilt(y2,h1b,h1a);
    end
    
elseif level == -1,
    
    for f = t3,
        % Do odd top-level filters on rows.
        y = colfilter(x(t1,s2,f).',g0o) + colfilter(x(t1,s2b,f).',g1o);
        % Do odd top-level filters on columns.
        x(s1,s2,f) = colfilter(y(:,s1).',g0o) + colfilter(y(:,s1b).',g1o);
    end
    
    for f = s2,
        % Do odd top-level filters on 3rd dim.
        y = reshape(x(s1,f,t3),sr(1),sx(3)).';
        x(s1,f,s3) = (colfilter(y(s3,:),g0o) + colfilter(y(s3b,:),g1o)).';
    end
    
elseif level <= -2,
    
    for f = t3,
        % Do even Qshift filters on rows.
        y = colifilt(x(t1,s2,f).',g0b,g0a) + colifilt(x(t1,s2b,f).',g1b,g1a);
        % Do even Qshift filters on columns.
        x(t1,t2,f) = colifilt(y(:,s1).',g0b,g0a) + colifilt(y(:,s1b).',g1b,g1a);
    end

    for f = t2,
        y = reshape(x(t1,f,t3),sx([1 3])).';
        % Do even Qshift filters on 3rd dim.
        x(t1,f,t3) = (colifilt(y(s3,:),g0b,g0a) + colifilt(y(s3b,:),g1b,g1a)).';
    end
    
else
   disp('Illegal level');
end

tk = toc;
fprintf(1,'toc = %.2f sec\n',tk);

return;

         
