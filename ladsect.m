function yy = ladsect(xx,h4)

% function yy = ladsect(xx,h4)
% Implement one section of ladder filtering on the 2 columns of xx.
% Columns 1 and 2 of xx are the z polynomials of signal x1 and x2.
% Columns 1 to 4 of h4 are the z polynomials of the 4 filters which
% convolve with x1 and x2 to give y1 and y2.

yy(:,1) = conv(xx(:,1),h4(:,1)) + conv(xx(:,2),h4(:,2));
yy(:,2) = conv(xx(:,1),h4(:,3)) + conv(xx(:,2),h4(:,4));



