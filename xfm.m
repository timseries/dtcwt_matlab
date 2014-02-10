function h = xfm(ht,m)
% function h = xfm(ht,m)
% Function to compute the 1-D filter coefs of the transformed filter h(z)
% from ht(Z) where Z = m(z).
% All polynomials are column vectors in powers of Z or z in descending order.
% m(z) is assumed to be of odd length with coef of z^0 as centre term.
M = length(m);
N = length(ht);
deg = N-1;
nul = zeros((M-1)/2,1);
h = ht(N);
Z = 1;
for i = 1:deg,
   Z =	conv(Z,m);
   h = [nul; h; nul] + ht(N-i) * Z;
end
