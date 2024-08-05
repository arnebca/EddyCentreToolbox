function y = summiss(x)
%SUMMISS Sum of elements with missing data.
%  Y = SUMMISS(X) returns the sum of each column of X as a row vector
%  where missing data values are encoded as NaNs. For vectors, SUMMISS(X)
%  returns the sum of the elements in X.

%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/09/27 14:45:28 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995

if ~isempty(x)
  [m,n] = size(x);
  y = zeros(m,n);
  valid = ~isnan(x);
  y(valid) = x(valid);
  y = sum(y);
else
  y=nan;
end

