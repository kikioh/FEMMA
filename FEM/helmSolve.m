% FEM solution to the helmholtz Filter equation. xf = R^2 * d^2(xf)/d^2(xc)
% + k_1 * (xc - xc0)
% Syntax:
%
%   [xf,grad] = helmholtzSolve(K,T,xc,k1,xc0)
%
%   Parameters
%       
%       R       - the stiffness matrix
%       T       - time duration
%       xc      - input waveform
%       k1      - scale coefficient, the default value is 1
%       xc0     - mean value of xc, default value is 0.5
%
%   Output
%
%       xf       - row vector, filtered waveform
%       grad     - N by N matrix, gradient of filtered waveform 
%                  to input waveform, i.e.,
%       [df1/dc1, df1/dc2,  df1/dc3,  .... df1/dcN;
%        df2/dc1, df2/dc2,  df2/dc3,  .... df2/dcN;
%        ...
%        dfN/dc1, dfN/dc2,  dfN/dc3,  .... dfN/dcN]
%
% mengjia.he@kit.edu, 2024.11.28

function [xf,grad] = helmSolve(K,T,xc,k1,xc0)

% scale factor for smoothed variable xf
if ~exist('k1','var'), k1 = 1; end

% Default mean value of xc
if ~exist('xc0','var'), xc0 = 0; end

% transpose for row vector
if (size(xc,2) ~= 1) && (size(xc,1) == 1)
    xc = transpose(xc);
end

N = numel(xc);
dt = T/N;

% constract load vector
f = dt/2 * ([xc(1);(xc(1:N-1)+xc(2:N));xc(N)] - xc0*[1;2*ones(N-1,1);1]);

if nargout == 1
	xf = linEqSolve(K,f);
	xf = k1 * transpose(xf);
	
elseif nargout == 2 
	
	[xf,grad] = linEqSolve(K,f);
	xf = k1 * transpose(xf);
	grad = k1 * grad * dt/2;
	
end

end


% Solve a linear system 
function [x,grad] = linEqSolve(K,f)

% FEM mesh
nn = length(f);     % number of dofs
ne = nn-1;          % number of elements

% solve for nodal dofs
x=K\f;
x=x(1:ne);
    
% gradient of varaibels to the controls
if nargout == 2
 
    df_dx = dfdphi(ne);
    grad = K \ df_dx;
    grad = grad(1:ne,:);

end

end