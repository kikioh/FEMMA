% calculate element stiffness matrix and its gradient for solving LvN equatoin. 
% Using the 2-points linear shape function.
% Syntax:
%
%       [ke,keGrad] = eleMatrix_linear(control)
%
% Parameters
%
%		control.pulse_dt		- row vector, pulse time slices 
%       control.offsets         - cell array, resonance offset in Hz
%		control.off_ops			- cell array, offset spin operators corresponding to the offsets
%       control.pwr_level       - cell array, B1 amplitude in rad/s
%       control.operators       - cell array, spin operators matrix, 
%       control.waveform        - matrix, represent non-unit control waveform
% 		control.pulse_form		- string, 'xy' for Cartesian coordinate components, 'phase' for phase components
%
% Output
%
%       ke          - element stiffness matrix
%       keGrad     	- gradient of element stiffness matrix wrt controls
%
% mengjia.he@kit.edu, 2024.11.12


function [ke,keGrad] = eleMatrix_linear(control)

% extract control parameters
numP = numel(control.pulse_dt);					% number of time slices
waveform = control.waveform;					% pulse waveform  
dt = control.pulse_dt(1);						% time slice duration
dofNode = size(control.operators{1},1);			% dofs of each node
L = control.drifts;								% drift Hamiltonian matrix

% coefficient for building element stiffness matrix
CE  = [-1/2  1/2; -1/2 1/2];
CL = 1i*dt * [2/6  1/6; 1/6  2/6];

% preallocate a same-size unit matrix
E = repmat(eye(dofNode), 1, 1, numP);

% preallocate Cartesian waveform
if strcmp(control.pulse_form,'phase')
    numSpin = size(waveform,1);        

elseif strcmp(control.pulse_form,'xy')
    numSpin = size(waveform,1)/2;    
end

% % preallocate element stiffness matrix
% E = repmat(eye(dofNode), 1, 1, numP);
% L = zeros(dofNode,dofNode,numP,'like', control.operators{2}); 

% element sttiffness matrix and its gradient
switch control.pulse_form

    case 'phase' 
		
		for k = 1:numSpin
			
			% transform the waveform into 3D matrix
			cos_wave = reshape(cos(waveform(k,:)),1,1,[]);
			sin_wave = reshape(sin(waveform(k,:)),1,1,[]);
			
			% addup Hamiltonian for selected spin 
			L = L + 2*pi*control.offsets{k} * control.off_ops{k} +...
					control.pwr_levels{k} * (control.operators{2*k-1} .* cos_wave + control.operators{2*k} .* sin_wave);
		end

        % build block-wise element stiffness matrix
		ke = zeros(2*dofNode, 2*dofNode, numP, 'like', control.operators{2});
        for i = 1:2
            rows = (i-1)*dofNode + (1:dofNode);
            for j = 1:2
                cols = (j-1)*dofNode + (1:dofNode);
                ke(rows,cols,:) = CE(i,j)*E + CL(i,j)*L;
            end
		end
 
        % gradient of element stiffness matrix
        if nargout == 2

            % preallocate gradient answer
            keGrad = cell(numSpin, 1);
            for k=1:numSpin
				
				% transform the waveform into diagonal sparse matrix
				cos_wave = control.pwr_levels{k} * spdiags(cos(waveform(k,:)'), 0, numP, numP);
				sin_wave = control.pwr_levels{k} * spdiags(sin(waveform(k,:)'), 0, numP, numP);
			
				% temporary spin oprtators
				opX =  kron(CL,control.operators{2*k-1});
				opY =  kron(CL,control.operators{2*k}); 
				
                % gradient of element stiffness matrix
				keGrad{k} = kron(-sin_wave, opX)  + kron(cos_wave, opY);

            end
			
			% combine gradient matrix
			keGrad = blkdiag(keGrad{:});
        end

    case 'xy'

        for k=1:numSpin

            % expansion waveform
            cos_wave = reshape(waveform(2*k-1,:), 1, 1, []);
            sin_wave = reshape(waveform(2*k,:), 1, 1, []);

            % Implicit expansion Hamiltonian
            L = L + 2*pi*control.offsets{k} * control.off_ops{k} + ...
                    control.pwr_levels{k} * (control.operators{2*k-1} .* cos_wave +...
											 control.operators{2*k} .* sin_wave);
        end

        % build block-wise element stiffness matrix
		ke = zeros(2*dofNode, 2*dofNode, numP, 'like', control.operators{2});
        for i = 1:2
            rows = (i-1)*dofNode + (1:dofNode);
            for j = 1:2
                cols = (j-1)*dofNode + (1:dofNode);
                ke(rows,cols,:) = CE(i,j)*E + CL(i,j)*L;
            end
		end

        % gradient of element stiffness matrix
        if nargout == 2

			keGrad = cell(numSpin*2, 1);
            for k=1:numSpin
			
				% temporary spin oprtators
				opX =  control.pwr_levels{k} * kron(CL,control.operators{2*k-1});
				opY =  control.pwr_levels{k} * kron(CL,control.operators{2*k}); 
				
				% gradient of element stiffness matrix for selected channel 
				keGrad{2*k-1} = kron(speye(numP), opX);
				keGrad{2*k} = kron(speye(numP), opY); 				
				
            end
            
			% combine gradient matrix
			keGrad = blkdiag(keGrad{:});
			
        end
end

end