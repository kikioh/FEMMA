% calculate element stiffness matrix and its gradient for solving LvN equatoin. 
% Using the quadratic shape function.
% Syntax:
%
%       [ke,keGrad] = eleMatrix_quad(control)
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
%       ke          - sparse matrix, element stiffness matrix
%       keGrad     	- sparse matrix, gradient of element stiffness matrix wrt controls
%
% mengjia.he@kit.edu, 2024.11.12


function [ke,keGrad] = eleMatrix_quad(control)

% number of outputs
nOutput = nargout;

% extract operators
numP = numel(control.pulse_dt);
waveform = control.waveform;
dt = control.pulse_dt(1)/2;
dofNode = size(control.operators{1},1);
dofElem = 3*dofNode;

% coefficient for building element stiffness matrix
CE  = [-1/2  2/3 -1/6;
	   -2/3  0    2/3;
		1/6 -2/3  1/2];

CL = 1i*dt* [ 4/15  2/15 -1/15;
			  2/15  16/15 2/15;
			 -1/15  2/15  4/15];  

% preallocate Cartesian waveform
if strcmp(control.pulse_form,'phase')
    numSpin = size(waveform,1);        

elseif strcmp(control.pulse_form,'xy')
    numSpin = size(waveform,1)/2;    
end

% extract drift Hamiltonian
L = control.drifts;

% preallocate a same-size unit matrix
E = repmat(eye(dofNode), 1, 1, numP);

% element sttiffness matrix and its gradient
switch control.pulse_form

    case 'phase'

		
        for k=1:numSpin

			% transform the waveform into 3D matrix
			cos_wave = reshape(cos(waveform(k,:)),1,1,[]);
			sin_wave = reshape(sin(waveform(k,:)),1,1,[]);

            % Implicit expansion Hamiltonian
            L = L + 2*pi*control.offsets{k} * control.off_ops{k} + ...
                    control.pwr_levels{k} * (control.operators{2*k-1} .* cos_wave +...
											control.operators{2*k} .* sin_wave);
        end
	
		% build element stiffness matrix 
		ke = zeros(dofElem,dofElem,numP,'like', control.operators{2});
		for i = 1:3
			rows = (i-1)*dofNode + (1:dofNode);
			for j = 1:3
				cols = (j-1)*dofNode + (1:dofNode);
				ke(rows,cols,:) = CE(i,j)*E + CL(i,j)*L;
			end
		end

        % gradient of element stiffness matrix
        if nOutput == 2
		
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

		% build element stiffness matrix 
		ke = zeros(dofElem,dofElem,numP,'like', control.operators{2});
		for i = 1:3
			rows = (i-1)*dofNode + (1:dofNode);
			for j = 1:3
				cols = (j-1)*dofNode + (1:dofNode);
				ke(rows,cols,:) = CE(i,j)*E + CL(i,j)*L;
			end
		end

        % gradient of element stiffness matrix
        if nOutput == 2

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