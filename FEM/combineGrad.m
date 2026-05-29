% Combine the gradient from LvN equation and the gradient from variables
% precessing function. Syntax:
%
%   drho_dx = combineGrad(drho_dxs, Dx)
%
% Parameters:
%
%   drho_dxs    - dofNode * nv matrix, gradient from LvN equation
%
%   Dx          - 1*numChn cell array, gradient from  variables
%               precessing equation
%
% Outputs:
%
%   drho_dx     - dofNode * nv matrix, gradient of rho to control variables
%
% mengjia.he@kit.edu, 2025.03.06

function dfdx = combineGrad(drho_dxs, Dx)

% extrac control parameters
numChn = numel(Dx);
[dofNode,nv] = size(drho_dxs);

% reshape gradient from LvN equation
drho_dxs = reshape(full(drho_dxs), dofNode, numChn, []);
drho_dxs = permute(drho_dxs, [1, 3, 2]);

% combine LvN gradient with variable process gradient, using chain rule
drho_dx = zeros(size(drho_dxs));
for k = 1 : numChn

    drho_dx(:,:,k) = drho_dxs(:,:,k) * Dx{k};

end

% reshape gradient of weight rho with respect to control variables
drho_dx = reshape(permute(drho_dx, [1, 3, 2]),dofNode,nv);

% gradient of fidelity with respect to control variables
% if dofNode > 1
	% dfdx = real(sum(drho_dx,1));
	
	% else
	dfdx = real(drho_dx);
% end	

end

