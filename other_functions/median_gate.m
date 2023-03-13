function [x,varargout] = median_gate(x,N,minpts,method,thrshld,varargin)
%MEDIAN_GATE Apply median gate to a vector
%   y = MEDIAN_GATE(x,N,minpts,method,K) replaces rejected elements of
%   x with NaN.  N array elements to either side of the test element
%   are used in the computation.  The gate is applied starting on
%   the center element of x and proceeds out to either side.
%
%   At least minpts must not have already been rejected, or the test
%   element is rejected as well.  Set minpts to 0 to turn this off.
%
%   Specify method as 'fixed', 'StdDev', or 'MADev.'  K specifies the
%   rejection threshold.  For 'fixed' a gain of K rejects elements
%   whose distance from the median is larger than K.  For 'StdDev' a
%   gain of K rejects elements whose distance from the median is
%   larger than K times the standard deviation of the 2N+1 sample set.
%   'MADev' works the same way except the deviation is computed using
%   mean absolute deviation instead of standard deviation.  This is a
%   more robust measure of sample deviation for samples containing
%   outliers.
%
%   [y,rjct_ptr] = MEDIAN_GATE( ... ) also returns an array of
%   indices to the rejected elements of x.
%
%   [ ... ] = MEDIAN_GATE( ... ,rjct_ptr) pre-rejects elements
%   specified in rjct_ptr.  This is useful in conjunction with
%   other gates or filters.  NaNs in x are treated as specifying a
%   rjct_ptr implicitly.
%
%   [ ... ] = MEDIAN_GATE( ... ,rjct_ptr,cx) specifies the first
%   element of x to apply the gate to.  Default is the central element.


% Apply input reject ptr
if nargin >= 6
  x(varargin{1}) = NaN;
end

% generate pattern of application
M = length(x);
if nargin >= 7
  cM = varargin{2};
else
  cM = ceil(M/2);
end
mm1 = cM:1:M;
mm2 = cM-1:-1:1;
mm(1:2:2*length(mm1)) = mm1;
mm(2:2:2*length(mm2)) = mm2;
mm = mm(find(mm)); % compress zeros

% apply median gate: output is x, not y
for m = mm
  y = sort(x(max(1,m-N):min(m+N,M)));
  y = y(find(~isnan(y)));
  npts = length(y);
  if npts < minpts | npts == 0
    x(m) = NaN;
  else
    y_med = y(ceil(npts/2));
    switch method
      case 'fixed'
	if abs(x(m)-y_med) > thrshld
	  x(m) = NaN;
	end
      case 'StdDev'
	if abs(x(m)-y_med) > thrshld*std(y)
	  x(m) = NaN;
	end
      case 'MADev'
	MADev = 1/(2*N+1)*sum(abs(y-y_med));
	if abs(x(m)-y_med) > thrshld*MADev
	  x(m) = NaN;
	end
      otherwise
	error('bad threshold option');   
    end
  end
end

% create output rjct_ptr
if nargout >= 2
  varargout{1} = find(isnan(x));
end
