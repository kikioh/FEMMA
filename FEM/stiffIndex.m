% Generate index matrix for assembling stiffness matrix. Syntax:
%
%   [rows,cols,rows_grad,cols_grad] = stiffIndex(nElements,nNodeDofs)
%
% Parameters:
%
%   nElement        - number of elements
%   nNodeDof        - dofs of each node
%
% Output
%
%   rows            -  row index matrix for global stiffness matrix
%   cols            -  column index matrix for global stiffness matrix
%   rows_grad       -  row index matrix for gradient stiffness matrix
%   cols_grad      	-  column index matrix for gradient stiffness matrix
%
% mengjia.he@kit.edu, 2025.01.21

function [rows,cols,rows_grad,cols_grad] = stiffIndex(nElems,dofNode,shapeOrder)

% Count the outputs
n_outputs=nargout;

% Default shape function
if ~exist('shapeOrder','var'), shapeOrder = 'linear'; end

% number of dof of element
switch shapeOrder

    case {'linear','hermite'}
        % element dofs for linear shape function
        dofElem = 2 * dofNode;
		
    case 'quadratic'
        dofElem = 3 * dofNode;
		
    case 'cubic'
        dofElem = 4 * dofNode;
end

% prelocate index matrix for stiffness matrix
rows = zeros(dofElem,dofElem,nElems);
cols = zeros(dofElem,dofElem,nElems);

if n_outputs <=2
    for m = 1:nElems
        % get the global scatter vector
%         sctr = dofNode*m-dofNode+1 : dofNode*m+dofNode;
        sctr =(dofElem-dofNode)*(m-1)+1:(dofElem-dofNode)*(m-1)+dofElem;
        rows(:,:,m) = sctr' * ones(1, dofElem);
        cols(:,:,m) = ones(dofElem, 1) * sctr;
    end
	
	rows = int32(rows);
	cols = int32(cols);

else

    % prelocate index matrix for gradient stiffness matrix
    rows_grad = zeros(dofElem*2,dofElem*2,nElems);
    cols_grad = zeros(dofElem*2,dofElem*2,nElems);

    for m = 1:nElems

        % get the global scatter vector
        sctr = dofNode*m - dofNode + 1 : dofNode*m + dofNode;
        rows(:,:,m) = sctr' * ones(1, dofElem);
        cols(:,:,m) = ones(dofElem, 1) * sctr;
        sctr_grad = dofElem*2*m-dofElem*2+1 : dofElem*2*m;
        rows_grad(:,:,m) = sctr_grad' * ones(1, dofElem*2);
        cols_grad(:,:,m)= ones(dofElem*2, 1) * sctr_grad;

    end
end

end