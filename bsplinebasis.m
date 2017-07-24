## B = bsplinebasis(degree,knots,ts)
##
## Evaluate basis for a spline of given degree and knots, at a set of parameter
## values ts.
##
## Each column corresponds to a control point and each row corresponds to one
## parameter value such that, for a column vector of control points c, the
## spline can be evaluated as:
##
## s = B.c
##
## Matthew Chapman <contact@zmatt.net>, based on algorithm in
##    C code for An Introduction to NURBS, David F. Rogers, Section 3.5

function B = bsplinebasis(degree,knots,ts)

order = degree + 1;
nknots = length(knots);
npoints = nknots-order;
B = [];

for i = [1:length(ts)]
  t = ts(i);
  if t == knots(nknots)
    # special case: if t is last knot, treat it as if it were
    # in the previous span, otherwise it would not be in any
    j = find(knots!=t,1,"last");
    searcht = knots(j);
  else
    searcht = t;
  endif

  # calculate 1st order basis functions
  # 1 if in knot span, 0 if not
  for j = [1:nknots-1]
    temp(j) = double(searcht >= knots(j) && searcht < knots(j+1));
  endfor

  for k = [2:order]
    # recursively calculate next order basis functions
    # by a linear combination of temp(j) and temp(j+1)
    for j = [1:nknots-k]
      d = e = 0;
      if temp(j) != 0
        d = ((t-knots(j))*temp(j))/(knots(j+k-1)-knots(j));
      endif
      if temp(j+1) != 0
        e = ((knots(j+k)-t)*temp(j+1))/(knots(j+k)-knots(j+1));
      endif
      temp(j) = d+e;
    endfor
  endfor

  B = [B; temp(1:npoints)];
endfor

endfunction
