function [h0,h1,g0,g1,H_lad,g_lad] = wavegen(wtype,param)

% function [h0,h1,g0,g1,H_lad,g_lad] = wavegen(wtype,param)
% Generate wavelet filters, where h0,h1 are low and high analysis filters
% (column vectors) and g0,g1 are low and high synthesis filters.
% The wavelet type is specified by the string wtype.
% param is a vector of parameters specific to each wavelet type.
% h0,g0 are adjusted for unit dc gain, and g1,h1 are the conjugate mirrors
% of these with unit gain at z = -1.
% For biorthogonal filters, g0 is the smoother filter.
% H_lad, g_lad (if required) are the H matrix and gain matrix for LADFIL.M .
% For odd-length filters, h1 introduces 1 sample less delay than h0.
% For even-length filters, h0 and h1 introduce the same delay.
%
% Currently supported values for wtype:
%   'Haar' = 2,2 tap Haar wavelet.
%   'LeGall' = 5,3 tap LeGall wavelet.
%   '6_10' = simple near symmetric even filters with binary coefs.
%   '6_10cplx' = version of 6,10 even filters for conversion from real to complex.
%   'Antonini' = 9,7 tap FBI filters.
%   'Daub44' = Daubechies 4,4-tap simplest wavelets
%   'Daub88' = Daubechies 8,8-tap near-symmetric wavelets.
%   'near-sym' = 5,7 and 13,19 tap near-symmetric filters:
%     param(1) = c (3.5), param(2:3) = [ma mb] ([3/16 0]).
%   'Qshift' = quarter sample shift orthogonal even-length filters:
%     param(1) = no of taps (must be even, typ 12 to 18),
%     param(2) = start of filter stop-band (typ 0.36).

if nargin < 2, param = 0; end

if strcmp(wtype,'Haar'),
% Start with the Haar wavelet.
  h0 = [1;1];
  g0 = h0;

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,  
% Generate ladder filter matrix and gain matrix.
    H_lad = [0.5; -1];
    g_lad = [1 0;0 -0.5]';
  end

elseif strcmp(wtype,'LeGall'),
% Now do the Le-Gall wavelet.
  h0 = [-1 2 6 2 -1]';
  g0 = [1 2 1]';

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,
% Generate ladder filter matrix and gain matrix.
    [H_lad,g_lad] = ladgen([h0 [h1;0;0]]);
  end

elseif strcmp(wtype,'6_10'),
% Calculate (6,10)-tap wavelet filters. 
% Z factors:
% ht=conv([1 3 3 1]',[2 -1]')
% conv(ht,[4 -2 -1 1]')
% h0=conv(xfm([-1 2]',[1 0 1]'/2),[1 3 3 1]')
% g0=conv(xfm([1 -1 -2 4]',[1 0 1]'/2),[1 3 3 1]')
  h0 = conv([1 3 3 1]',[-1 4 -1]');
  g0 = conv([1 3 3 1]',[1 -2 -5 28 -5 -2 1]');

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,  
% Generate ladder filter matrix and gain matrix.
% This gives 10-tap h0e and 6-tap h1e.
% He_lad = [0 -0.5 0; 0 1 0; [1 0 -1; 1 0 -1]/8];
% ge_lad = [1 0.5]';
% This gives 6-tap h0e and 10-tap h1e (the 10-tap g0e is smoother).
%    H_lad = [0 0.5 0; 0 -1 0; 0 0 0; [-1 0 1; -1 0 1]/8];
%    g_lad = [1 0;0 -0.5]';

    [H_lad,g_lad] = ladgen(beside(h0,h1));
  end

elseif strcmp(wtype,'6_10cplx'),
% Calculate (6,10)-tap wavelet filters for conversion from real to complex. 
  h0 = [1 -1 8 8 -1 1].';
  g0 = [1 1 8 -8 62 62 -8 8 1 1].';

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,  
% Generate ladder filter matrix and gain matrix.
% This gives 6-tap h0e and 10-tap h1e (the 10-tap g0e is smoother).
%    H_lad = [0 0.5 0; 0 -1 0; 0 0 0; [1 0 -1; 1 0 -1]/8];
%    g_lad = [1 -0.5]';
    [H_lad,g_lad] = ladgen(beside(h0,h1));
  end

elseif strcmp(wtype,'Antonini'),
% Antonini 9,7-tap wavelets (FBI standard).
% This polynomial makes even terms of P(z) zero when other factors 
% are (1+z)^8.
  r2 = roots([5 -40 131 -208 131 -40 5]);

% 4 roots at z=-1 and the 4 complex roots from r2.
  h0 = real(poly([-1; -1; -1; -1; r2(2:5)])');
% 4 roots at z=-1 and the 2 real roots from r2.
  g0 = real(poly([-1; -1; -1; -1; r2([1 6])])');

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,
% Generate ladder filter matrix and gain matrix.
    [H_lad,g_lad] = ladgen([h0 [h1;0;0]]);
  end

elseif strcmp(wtype,'Daub44'),
% Daubechies 4,4-tap simplest wavelets.
% This polynomial makes even terms of P(z) zero when other factors 
% are (1+z)^4.
  r2 = roots([-1 4 -1]);

% 4 roots at z=-1 and the 4 complex roots from r2.
  h0 = real(poly([-1; -1; r2(1)])');
% 4 roots at z=-1 and the 2 real roots from r2.
  g0 = real(poly([-1; -1; r2(2)])');

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,
% Generate ladder filter matrix and gain matrix.
    [H_lad,g_lad] = ladgen([h0 h1]);
  end

elseif strcmp(wtype,'Daub66'),
% Daubechies 4,4-tap simplest wavelets.
% This polynomial makes even terms of P(z) zero when other factors 
% are (1+z)^4.
  r2 = roots([3 -18 38 -18 3]);

% 4 roots at z=-1 and the 4 complex roots from r2.
  h0 = real(poly([-1; -1; -1; r2([1 2])])');
% 4 roots at z=-1 and the 2 real roots from r2.
  g0 = real(poly([-1; -1; -1; r2([3 4])])');

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,
% Generate ladder filter matrix and gain matrix.
    [H_lad,g_lad] = ladgen([h0 h1]);
  end

elseif strcmp(wtype,'Daub88'),
% Daubechies 8,8-tap near-symmetric wavelets.
% This polynomial makes even terms of P(z) zero when other factors 
% are (1+z)^8.
  r2 = roots([5 -40 131 -208 131 -40 5]);

% 4 roots at z=-1 and the 4 complex roots from r2.
  h0 = real(poly([-1; -1; -1; -1; r2([1 4 5])])');
% 4 roots at z=-1 and the 2 real roots from r2.
  g0 = real(poly([-1; -1; -1; -1; r2([2 3 6])])');

  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,
% Generate ladder filter matrix and gain matrix.
    [H_lad,g_lad] = ladgen([h0 h1]);
  end

elseif strcmp(wtype,'near-sym'),
% Transformation of Variables odd-tap filters:
% Now do the 5,7-tap and 13,19-tap Transformation of variables filters.
  if nargin < 2, param = [3.5 0 0]; end

  c = -param(1); 
  a = 2*c + 2/(2 + c);
  b = -(2 + c);
  Ht = conv([1 1]',-[1 b a]');
  Ft = conv([1 1]',-[1 c]');

% Define the transformation function Z = M(z) (even coefs must be 0).
% for the 13,19-tap filters.
  ma = param(2);  % 3/16 gives a good smooth filter.
  mb = param(3);  % 0 for 4-tap M.
  if ma > 0,
    if mb > 0,
      M = [mb 0 -ma 0 1+ma-mb 0 1+ma-mb 0 -ma 0 mb]'/2;
    else M = [-ma 0 1+ma 0 1+ma 0 -ma]'/2;
    end
  else M = [1 0 1]'/2;
  end

% Calculate the transformed low and high synthesis filters, lo(z) and hi(z).
  h0 = xfm(Ft,M);
  g0 = xfm(Ht,M);
  h0 = h0 / sum(h0);
  g0 = g0 / sum(g0);
  h1 = g0 .* cumprod(-ones(size(g0)));
  g1 = -h0 .* cumprod(-ones(size(h0)));

  if nargout > 4,
% Generate ladder filter matrix and gain matrix.
    [H_lad,g_lad] = ladgen([[0;0;0;0;h0;0;0] h1]);
  end
  
elseif strcmp(wtype,'Qshift') & (param(1) > 8),
% Quarter-sample shift CWT design based on a linear-phase
% 2x oversampled even-length filter (length 2N).
% This version minimises the energy over a defined frequency range.
% It works with the minimum number of unknowns and specifies
% 2 predefined zeros at z=-1.
  if nargin < 2, param = [14 0.36]; end

  N = param(1); % No of filter taps (must be even, typ 12 to 18).
  fmin = param(2); % Bottom of stopband (typ 0.36).
  if rem(N,2)>0, 
     error('Wavegen: Qshift filter order must be even!');
  end

% Design initial guess at prototype h filter.
  np = 16;
  mp = 4; 
  x = [0:np]' * (1/np);
  theta = [x.^mp; 2-x(np:-1:2).^mp] * (pi/4);
  Hp = [cos(theta); zeros(2*np,1)];  % Prototype freq response.

% Apply half-sample phase shift and do IFFT.
  Hp = Hp .* exp([0:(4*np-1)]' * (j*pi/(8*np)));
  hp = real(ifft([Hp; 0; conj(Hp((4*np):-1:2))]));
% Truncate impulse reponse to 2N taps and ensure it is symmetric.
  h = 2 * hp([N:-1:1  1:N]);
%  H0 = fft(h,8*np);

% Test PR accuracy of prototype.
%  pf=conv(h,h).'*100;
%  mat2scr([0 pf(4:4:4*(N-1))],'%9.4f');

% Loop to find h which satisfies PR
% and minimises energy of the frequency components defined by f0.
  mag = 20;
  fmag = 20;
  f0 = pi*(fmin + (1-fmin)*[0:0.5:N].'/N);
  fzero = exp(j*f0*([1:length(h)] - 1));
  fzero = [real(fzero);  imag(fzero)] * fmag;
  for loop=1:20,
   oldh = h;
   % fcorr corrects for the process    h = (h + oldh)/2   when minimising the
   % frequency terms.
   fcorr = fzero * oldh;
   cmat = toeplitz([h; zeros(2*N-1,1)],[h(1) zeros(1,2*N-1)]); % Convolution matrix.
   cm2 = [cmat([4:4:2*N],:)*mag; fzero]; % PR terms from conv. matrix + freq response rows.
   t = 1:(size(cm2,2)-4);
   cm2 = cm2(:,t) + 2*cm2(:,t+2) + cm2(:,t+4); % Convolve cm2 with [1 0 2 0 1]
   cm2 = cm2(:,1:(N-2)) + cm2(:,((N-2):-1:1)+N-2); % Combine columns for symmetric h coefs.
   h2 = cm2\[zeros(N/2-1,1); mag; -fcorr]; % Find h2 for best fit to PR and freq responses.
   h = h2([1:(N-2)  (N-2):-1:1]);         % Form symmetric h from h2.
   h = conv(h,[1 0 2 0 1]'); % Convolve h with [1 0 2 0 1].
   h = (h + oldh)/2;         % Average the old and new h's.
% The following line can force PR at each iteration but this slows convergence.
%   hpr = prforce(h(1:2:2*N)); h = reshape([hpr hpr(N:-1:1)].',2*N,1);
%   pf=conv(h,h).'*100;       % Test PR for new h.
%   mat2scr([mag/1000 sum(abs(h-oldh))*100 pf(4:4:4*(N-1))],'%9.4f');
   mag = mag*2;              % Increase mag to get progressively better PR fit.
  end
  
  h0 = h(2:2:2*N);
  h0 = prforce(h0);  % Force PR on final result.
  h0 = h0 / sum(h0);
  
  g0 = h0(N:-1:1);   % g0 is time reverse of h0.
  h1 = -g0 .* cumprod(-ones(size(g0)));  % h1 and g1 are quadrature mirrors
  g1 = h0 .* cumprod(-ones(size(h0))); % of g0 and h0 respectively.

  if nargout > 4,
% Generate ladder filter matrix and gain matrix.
    [H_lad,g_lad] = ladgen([h0 h1],3e-6);
  end
  
elseif strcmp(wtype,'Qshift'),
   % Theta(1) = pi/4 to give one zero at z=-1.
   % Theta(2) = Theta(1) to give 2 zero terms at h([2 7]) ie 6-tap complexity.
   % Theta(3:4) adjusted for best wavelet shape and smoothness.
   Theta = [1 1 -0.81 -1.62]*pi/4;
   % Alternative Theta's:
   % Theta = [1 1.07 -0.82 -1.64]*pi/4; % for lowest sidelobes in fig 2.
	% Theta = [1 1.03 -0.77 -1.62]*pi/4; % for smoothest wavelets at levels 3 on.

   % Calculate the impulse responses of the lattice network with cumulative 
	% rotations to the output of Theta(:).
	theta = -diff([Theta 0]);
	xx = [1 0;0 1];
	i = 1;
	c = cos(theta(i)); s = sin(theta(i));
	xx = ladsect(xx,[c s -s c]);  % Basic rotation.
	for i = 2:length(theta);
		c = cos(theta(i)); s = sin(theta(i));
		% Rotation with delay by 2 samples in lower input.
   	xx = ladsect(xx,[c 0 -s 0; 0 0 0 0; 0 s 0 c]); 
   end
   
   h0 = [xx(:,1);0;0] / sum(xx(:,1));
   N = length(h0);  
	g0 = h0(N:-1:1);   % g0 is time reverse of h0.
	h1 = -g0 .* cumprod(-ones(size(g0)));  % h1 and g1 are quadrature mirrors
	g1 = h0 .* cumprod(-ones(size(h0))); % of g0 and h0 respectively.

	if nargout > 4,
	% Generate ladder filter matrix and gain matrix.
   	[H_lad,g_lad] = ladgen([h0 h1],3e-6);
	end
  
else
  disp('Supported wtypes are: Haar, LeGall, 6_10, Antonini, near-sym, Qshift.');
  error('WAVEGEN.M: Unsupported wtype!');
end

return;




