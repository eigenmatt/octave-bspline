% [cpts,knots,degree] = bsplinederiv(degree,cpts,knots)
%
% Given a spline with given degree, control points and knots, return the
% spline corresponding to its derivative.  Note that the degree is always
% one less than the input degree.
%
% Matthew Chapman <contact@zmatt.net>, based on algorithm in
%   https://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-derv.html

function [cpts,knots,degree] = bsplinederiv(degree,cpts,knots)

% calculate new control points
for j = 1:length(cpts)-1
  alpha = degree/(knots(j+degree+1)-knots(j+1));
  cpts(j) = alpha * (cpts(j+1)-cpts(j));
end

cpts(end) = [];
knots = knots(2:end-1);
degree = degree-1;

end
