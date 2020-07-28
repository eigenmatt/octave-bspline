% vals = bsplineeval(degree,cpts,knots,ts)
%
% Given a spline of given degree, control points and knots, evaluate it at
% a set of parameter values ts.
%
% This is just a convenience wrapper around the bsplinebasis function.
%
% Matthew Chapman <contact@zmatt.net>

function vals = bsplineeval(degree,cpts,knots,ts)

vals = bsplinebasis(degree,knots,ts)*cpts;

end