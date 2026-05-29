% prepare parameters for FEM. Syntax:
%
%    fem = meshFEM(control,equationType,shapeOrder,R)
%
% Parameters:
%
%   control                 - control structure specifying the pulse 
%
%   equationType            - 'LvN' or 'Helmholtz'
%
%	shapeOrder (optional)	- name of FEM element, can be 'linear', 'quadratic' or 'hermite'
%
%	R						- Radius for Helmholtz filter
%
% Outputs:
%
%   fem.nNodes              - number of nodes
%   fem.nElems              - number of elements
%   fem.T                   - total time length
%   fem.dofNode             - number of dofs of each node
%	fem.K					- stiffness matrix for Helmholtz filter
%	fem.KIndex				- index of stiffness matrix for LvN equation
%
% mengjia.he@kit.edu, 2025.01.14

function [fem,control] = meshFEM(control,equationType,shapeOrder,R)

% extract time slice
T = sum(control.pulse_dt);
nElems = numel(control.pulse_dt);

% Default shape function
if ~exist('shapeOrder','var'), shapeOrder = 'linear'; end

% Default waveform type, piecewise-constant
if ~isfield(control,'integrator'), control.integrator = 'rectangle'; end

% Default radius for Helmholtz filter
if strcmp(equationType, 'Helmholtz') && (~exist('R','var')), R = T/130; end

% Only linear functions for bloch equation
if strcmp(equationType, 'bloch') && (~strcmp(shapeOrder,'linear'))
	error('Bloch equation only supports linear functions');
end

% number of dofs of each node
switch  equationType

    case 'Helmholtz'
		fem.R = R;
	
		switch shapeOrder
		
			case 'linear'  
			
				nNodes = nElems+1;  
				dofNode = 1;
			
			otherwise 
			
				error('Helmholtz filter only support linear shape functions');
				
		end

        % element stiffness matrix, independent on coordinates
		dt = T/nElems;
        S_elem = dt * [1/3, 1/6; 1/6, 1/3] + R^2/dt * [1,-1;-1,1];
		
		switch control.integrator
			
			case 'rectangle'
			
			    % combine element stiffness matrix
				S = repmat(S_elem, 1, 1, nElems);

				% assemble global stiffness matrix
				[rows,cols] = stiffIndex(nElems,dofNode);
				fem.K = sparse(rows(:), cols(:), S(:), nNodes*dofNode, nNodes*dofNode);
			
			case 'trapezium'
			
				% combine element stiffness matrix
				S = repmat(S_elem, 1, 1, nElems+1);

				% assemble global stiffness matrix
				[rows,cols] = stiffIndex(nElems+1,dofNode);
				fem.K = sparse(rows(:), cols(:), S(:), (1+nNodes)*dofNode, (1+nNodes)*dofNode);
		
		end

    case 'LvN'
		
		switch shapeOrder

			case 'linear'

				nNodes = nElems+1;  
				dofNode = size(control.operators{1},1);
				% dofElem = 2*dofNode;
				
				% gererate nonzero elements for matlab nonzeros function
				% control.ke_reg = 1e-10 * kron(speye(nElems), kron(ones(2), ones(dofNode)));
			
			case 'quadratic'
			
				nNodes = 2*nElems+1;  
				dofNode = size(control.operators{1},1);
				% dofElem = 3*dofNode;
				
				% gererate nonzero elements for matlab nonzeros function
				% control.ke_reg = 1e-10 * kron(speye(nElems), kron(ones(3), ones(dofNode)));
			
			case 'cubic'
			
				nNodes = 3*nElems+1;
				dofNode = size(control.operators{1},1);
				% dofElem = 4*dofNode;
				
				% gererate nonzero elements for matlab nonzeros function
				% control.ke_reg = 1e-10 * kron(speye(nElems), kron(ones(4), ones(dofNode)));
			
			case 'hermite'
			
				nNodes = nElems+1; 
				dofNode = 2*size(control.operators{1},1);
				% dofElem = 2*dofNode;
				
				% gererate nonzero elements for matlab nonzeros function
				% control.ke_reg = 1e-10 * kron(speye(nElems), kron(ones(4), ones(size(control.operators{1}))));

		end

        % index matrix for global stiffness matrix
        [fem.KIndex.rows,fem.KIndex.cols] = stiffIndex(nElems,dofNode,shapeOrder);  
%       fem.KIndex.rows_grad = rows_grad; fem.KIndex.cols_grad = cols_grad;

		% index matrix for assemble global force
		% nc = size(control.waveform,1);
		% control.dofs_matrix = dofIndex(dofNode, dofElem, nElems, nc, 'strided');
		
	case 'bloch'
	
		nNodes = nElems+1;  
		dofNode = 3;
		
		% index matrix for global stiffness matrix
        [fem.KIndex.rows,fem.KIndex.cols] = stiffIndex(nElems,dofNode,'linear');		 	
		
end

% assemble results to a structure
fem.nElems = nElems;
fem.T = T;
fem.equation = equationType;
fem.shapeOrder = shapeOrder;
fem.nNodes = nNodes;
fem.dofNode = dofNode;

end

