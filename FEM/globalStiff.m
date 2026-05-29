% assembly the global stiffness matrix and its gradient. Syntax:
%
%       [K,f,Kgrad,Kgrad_R] = globalStiff(fem,control)
%
% Parameters:
%
%		fem
%       	fem.T                  	- total time
%			fem.nEmems				- number of elements		
%       	fem.nNodes           	- number of nodes	
%       	fem.dofNode          	- dofs of node
%			fem.KIndex (optional)	- Index matrix for assembly K, pre-computed for 'LvN' equation
%
%		control
%       	control.offsets         - 1*numSpin cell array, resonance offset in Hz
%       	control.pwr_level      	- 1*numSpin cell array, B1 amplitude for all control channels in rad/s
%       	control.waveform      	- numChn by numP matrix, represent phase or xy waveform
%
% Output:
%
%       K                   	- global stiffness matrix
%       Kgrad               	- global stiffness matrix gradient wrt controls
%       Kgrad_R (optional)      - right-side gradient of global stiffness matrix wrt controls if Hemerite functions is used
%
% mengjia.he@kit.edu, 2025.09.04

function [K,f,Kgrad,Kgrad_R] = globalStiff(fem,control)

% uniform parameters
dt = fem.T/fem.nElems;              % slice duratoin
nElems = fem.nElems;                % number of elements
nNodes = fem.nNodes;                % number of nodes
dofNode = fem.dofNode;              % dofs of each node
dof = nNodes * dofNode;				% total dofs

switch  fem.equation

	case 'LvN'
		
		% all-zero load vector
		f = sparse(dof,1);
		
		% preallocate global stiffness matrix
		if nargout == 2

			% element stiffness matrix
			switch fem.shapeOrder

				case 'linear', ke = eleMatrix_linear(control);
				
				case 'quadratic', ke = eleMatrix_quad(control);
				
				case 'hermite', ke = eleMatrix_hermite(control);
				
			end

		elseif nargout >= 3

			% element stiffness matrix and gradient
			switch fem.shapeOrder

				case 'linear', [ke,Kgrad] = eleMatrix_linear(control);
				
				case 'quadratic', [ke,Kgrad] = eleMatrix_quad(control); 
				
				case 'hermite', [ke,Kgrad,Kgrad_R] = eleMatrix_hermite(control); 
				
			end

		end
		
	case 'bloch'
	
		% preallocate global stiffness matrix
		if nargout == 2
		
			[ke,f] = eleMatrix_bloch(fem, control);
			
		elseif nargout == 3
			
			[ke,f,Kgrad] = eleMatrix_bloch(fem, control);
		
		else 
			error('piecewise-linear approximation is not supported for Bloch equation, use pwc instead');
		
		end

end

% assemble global stiffness matrix
rows = fem.KIndex.rows; cols = fem.KIndex.cols; 
K = sparse(rows(:), cols(:), ke(:), dofNode*nNodes, dofNode*nNodes);

end
