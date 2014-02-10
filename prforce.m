function hpr = prforce(h)

% function hpr = prforce(h)
% Solve for symmetric ladder filter to implement the filter h
% and its time-reverse quadrature mirror h1.
% Return the value of the perfect reconstruction version of h.
% h must be a column vector.
%
% Nick Kingsbury, Cambridge University, Feb 99.

% [h0a,h1a,g0a,g1a] = wavegen('Qshift',[14 0.36]);

n = length(h);
n2 = fix(n/2) - 1;

h1 = h(n:-1:1) .* cumprod(-ones(n,1));
hh=[h h1];
a = []; 
d = [];
err = [];

for i = 1:n2,
	sh = size(hh);
	t = 3:sh(1);
	a(i) = sum(prod(hh(1:2,:).')) / sum(hh(1:2,2).^2);
   if abs(a(i)) < 1,
      hh1 = hh * [1 a(i); -a(i) 1];
      hh = [hh1(t,1)  hh1(t-2,2)];
      d(i) = 1;
   else
      a(i) = 1/a(i);
      hh1 = hh * [1 -a(i); a(i) 1];
      hh = [hh1(t-2,1)  hh1(t,2)];
      d(i) = 2;
   end
   err(:,i) = hh1(1:2,d(i));
end

% hh2 = hh
% a,d,err

z2 = [0;0];
for i = n2:-1:1,
	sh = size(hh);
	t = 3:sh(1);
   if d(i) == 1,
      hh1 = [[z2;hh(:,1)]  [hh(:,2);z2]];
      hh = hh1 * [1 -a(i); a(i) 1];
   else
      hh1 = [[hh(:,1);z2]  [z2;hh(:,2)]];
      hh = hh1* [1 a(i); -a(i) 1];
   end
end

hpr = hh(:,1) / prod(1+a.*a);
return

