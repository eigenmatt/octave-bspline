% [cpts,knots] = bsplineinsert(degree,cpts,knots,t)
%
% Given a spline with given degree, control points and knots, insert a new
% knot at t and return new control points and knots.
%
% Matthew Chapman <contact@zmatt.net>, based on algorithm in
%   https://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/NURBS-knot-insert.html

function [cpts,knots] = bsplineinsert(degree,cpts,knots,t)

% find knot span that contains t
for k = 1:length(knots)-1
  if t >= knots(k) && t < knots(k+1)
    break
  end
end

% calculate new control points
newcpts = [];
for j = k-degree+1:k
  alpha = (t-knots(j))/(knots(j+degree)-knots(j));
  newcpts = [newcpts; (1-alpha)*cpts(j-1)+alpha*cpts(j)];
end

cpts = [cpts(1:k-degree); newcpts; cpts(k:end)];
knots = [knots(1:k) t knots(k+1:end)];

end
