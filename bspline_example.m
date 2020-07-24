%% B-spline example
%
% This is a simple example illustrating the creation of a B-spline basis,
% the definition of a particular spline as a linear combination of the
% basis elements, and the computation of its derivatives.

%% main parameters for example
degree = 3; % degree of the spline basis
knots = [0,0,0,0,0.25,0.5,0.75,1,1,1,1]; % knot vector
xrange = -0.3:0.01:1.3; % range of values over which we will evaluate the splines

%% define basis functions and plot them
basis = bsplinebasis(degree, knots, xrange);

clf;
subplot(1,2,1);
plot(xrange, basis)
title("B-spline basis")

%% plot a random example linear combination of the basis elements and its first two derivatives
cpts = rand(size(basis,2),1);
[d1cpts, d1knots, d1degree] = bsplinederiv(degree, cpts, knots);
[d2cpts, d2knots, d2degree] = bsplinederiv(d1degree, d1cpts, d1knots);

subplot(3,2,2);
plot(xrange, bsplineeval(degree, cpts, knots, xrange));
title("Random spline")
subplot(3,2,4);
plot(xrange, bsplineeval(d1degree, d1cpts, d1knots, xrange));
title("First derivative")
subplot(3,2,6);
plot(xrange, bsplineeval(d2degree, d2cpts, d2knots, xrange));
title("Second derivative")
