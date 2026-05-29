% Consistency enforcement
function grumbleMMA(m,n,xval,xmin,xmax,xold1,xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d)

%  xval  = Column vector with the current values of the variables x_j.
if ~isequal(size(xval), [n, 1])
    error('Current values of the variables must be a n by 1 column vector');
end

%  xmin  = Column vector with the lower bounds for the variables x_j.
if ~isequal(size(xmin), [n, 1])
    error('Lower bounds of the variables must be a n by 1 column vector');
end

%  xmax  = Column vector with the upper bounds for the variables x_j.
if ~isequal(size(xmax), [n, 1])
    error('Upper bounds of the variables must be a n by 1 column vector');
end

%  xold1 = xval, one iteration ago (provided that iter>1).
if ~isequal(size(xold1), [n, 1])
    error('xval of one iteration ago must be a n by 1 column vector');
end

%  xold2 = xval, two iterations ago (provided that iter>2).
if ~isequal(size(xold2), [n, 1])
	error('xval of 2 iteration ago must be a n by 1 column vector');
end

%  f0val = The value of the objective function f_0 at xval.
if ~isnumeric(f0val)
   error('f0val should be a number');
end

%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
if ~isequal(size(df0dx), [n, 1])
	error('df0dx must be a n by 1 column vector');
end

% df0dx2 = Column vector with the non-mixed second derivatives of the 
%          objective function f_0 with respect to the variables x_j,
%          calculated at xval. df0dx2(j) = the second derivative
%          of f_0 with respect to x_j (twice).
%          Important note: If second derivatives are not available,
%          simply let df0dx2 = 0*df0dx.
if ~isequal(size(df0dx2), [n, 1])
	error('df0dx2 must be a n by 1 column vector');
end

%  fval  = Column vector with the values of the constraint functions f_i,
%          calculated at xval.
if ~isequal(size(fval), [m, 1])
	error('df0dx2 must be a n by 1 column vector');
end

%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
if ~isequal(size(dfdx), [m, n])
	error('dfdx must be a n by 1 column vector');
end

%  dfdx2 = (m x n)-matrix with the non-mixed second derivatives of the
%          constraint functions f_i with respect to the variables x_j,
%          calculated at xval. dfdx2(i,j) = the second derivative
%          of f_i with respect to x_j (twice).
%          Important note: If second derivatives are not available,
%          simply let dfdx2 = 0*dfdx.
if ~isequal(size(dfdx2), [m, n])
	error('dfdx2 must be a n by 1 column vector');
end

%  low   = Column vector with the lower asymptotes from the previous
%          iteration (provided that iter>1).
if ~isequal(size(low), [n, 1])
	error('the lower asymptotes must be a n by 1 column vector');
end

%  upp   = Column vector with the upper asymptotes from the previous
%          iteration (provided that iter>1).
if ~isequal(size(upp), [n, 1])
	error('the upper asymptotes must be a n by 1 column vector');
end

%  a0    = The constants a_0 in the term a_0*z.
if ~isnumeric(a0)
   error('a0 should be a number');
end

%  a     = Column vector with the constants a_i in the terms a_i*z.
if ~isequal(size(a), [m, 1])
	error('constants a_i must be a m by 1 column vector');
end

%  c     = Column vector with the constants c_i in the terms c_i*y_i.
if ~isequal(size(c), [m, 1])
	error('constants c must be a m by 1 column vector');
end

%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
if ~isequal(size(d), [m, 1])
	error('constants d must be a m by 1 column vector');
end

end