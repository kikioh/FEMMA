% calculate element stiffness matrix and its gradient for solving LvN equatoin. 
% Using cubic Hermite shape function.
%
% Note: the contorls are defined on the nodes and the waveform is piecewise-linear, which is 
% refered to 2023_JMR_Simulation and design of shaped pulses beyond the piecewise-constant approximation.
%
% Syntax:
%
%       [ke,keGrad_L,keGrad_R] = eleMatrix_hermite(control)
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
% Output:
%
%       ke          	- element stiffness matrix
%       keGrad_L     	- gradient of left-side element stiffness matrix wrt controls
%       keGrad_R     	- gradient of right-side element stiffness matrix wrt controls
%
% mengjia.he@kit.edu, 2025.09.02


function [ke,keGrad_L,keGrad_R] = eleMatrix_hermite(control)

% number of outputs
nOutput = nargout;

% extract parameters
waveform = control.waveform;			% waveform in phase or xy
dt = control.pulse_dt(1);				% constant slice length
numP = numel(control.pulse_dt);			% number of time slices
numNode = size(waveform,2);				% number of nodes
nRho = size(control.operators{1},1);	% demension of spin vector 
dofNode = 2 * nRho;						% dof of node
dofElem = 2*dofNode;					% dof of element

% preallocate Cartesian waveform
if strcmp(control.pulse_form,'phase')
    numSpin = size(waveform,1);        

elseif strcmp(control.pulse_form,'xy')
    numSpin = size(waveform,1)/2;  
	numChn = 2*numSpin;
end

% coefficient for building element stiffness matrix
CE  = [-1/2  1/5   1/2  -1/5;
	   -1/5   0    1/5  -1/15;
	   -1/2 -1/5   1/2   1/5;
		1/5  1/15 -1/5   0];

CL1 = 1i*dt* [ 2/7   1/14   9/140  -1/30;
			   1/14  1/42   1/35   -1/70;
			   9/140 1/35   3/35   -1/30;
			   -1/30 -1/70  -1/30    1/70];

CL2 = 1i*dt* [ 3/35   1/30   9/140  -1/35;
			   1/30   1/70   1/30   -1/70;
			   9/140  1/30   2/7    -1/14;
			   -1/35  -1/70  -1/14    1/42];

% % preallocate Hamiltonian
% E = repmat(eye(nRho), 1, 1, numP);
% L1   = zeros(nRho,nRho,numP, 'like', control.operators{2}); 
% L2   = zeros(nRho,nRho,numP, 'like', control.operators{2});  

% extract drift Hamiltonian
L1 = control.drifts(:,:,1:end-1);
L2 = control.drifts(:,:,2:end);

% preallocate a same-size unit matrix
E = repmat(eye(nRho), 1, 1, numP);

% element sttiffness matrix and its gradient
switch control.pulse_form

    case 'phase'    

		% build L1 and L2
		for k = 1:numSpin
			
			% transform waveform to Cartesian coordinates
            cos_wave = reshape(cos(waveform(k,:)), 1, 1, []);
            sin_wave = reshape(sin(waveform(k,:)), 1, 1, []);
			
			% Implicit expansion Hamiltonian
            L1 = L1 + 2*pi*control.offsets{k} * control.off_ops{k} + ...
                    control.pwr_levels{k} * (control.operators{2*k-1} .* cos_wave(1:end-1) +...
                                            control.operators{2*k} .* sin_wave(1:end-1));
					
			L2 = L2 + 2*pi*control.offsets{k} * control.off_ops{k} + ...
                    control.pwr_levels{k} * (control.operators{2*k-1} .* cos_wave(2:end) +...
											control.operators{2*k} .* sin_wave(2:end));
		end

		% build element stiffness matrix 
		ke = zeros(dofElem,dofElem,numP,'like', control.operators{2});
		for i = 1:4
			rows = (i-1)*nRho + (1:nRho);
			for j = 1:4
				cols = (j-1)*nRho + (1:nRho);
				ke(rows,cols,:) = CE(i,j)*E + CL1(i,j)*L1 + CL2(i,j)*L2;
			end
		end

		% compute gradient if required
		if nOutput == 3
			
			% preallocate gradient answer
			keGrad_R = cell(numSpin, 1);
            keGrad_L = cell(numSpin, 1);
			
			for k = 1:numSpin
				
				% transform the waveform into diagonal sparse matrix
				cos_wave = control.pwr_levels{k} * spdiags(cos(waveform(k,:)'), 0, numP+1, numP+1);
				sin_wave = control.pwr_levels{k} * spdiags(sin(waveform(k,:)'), 0, numP+1, numP+1);
			
				% extend spin operators
				opX1 =  kron(CL1,control.operators{2*k-1}); opX2 =  kron(CL2,control.operators{2*k-1});
				opY1 =  kron(CL1,control.operators{2*k});   opY2 =  kron(CL2,control.operators{2*k}); 
				
                % gradient of element stiffness matrix
				keGrad_R{k} = kron(-sin_wave(1:end-1,1:end-1), opX1)  + kron(cos_wave(1:end-1,1:end-1), opY1);
				keGrad_L{k} = kron(-sin_wave(2:end,2:end), opX2)  + kron(cos_wave(2:end,2:end), opY2);
				
			end
			
			% combine gradient matrix
			keGrad_L = blkdiag(keGrad_L{:});
			keGrad_R = blkdiag(keGrad_R{:}); 
		end

    case 'xy'
		
        for k=1:numSpin

            % expansion waveform
            cos_wave = reshape(waveform(2*k-1,:), 1, 1, []);
            sin_wave = reshape(waveform(2*k,:), 1, 1, []);

            % Implicit expansion Hamiltonian
            L1 = L1 + 2*pi*control.offsets{k} * control.off_ops{k} + ...
                    control.pwr_levels{k} * (control.operators{2*k-1} .* cos_wave(1:end-1) +...
											 control.operators{2*k} .* sin_wave(1:end-1));
					
			L2 = L2 + 2*pi*control.offsets{k} * control.off_ops{k} + ...
                    control.pwr_levels{k} * (control.operators{2*k-1} .* cos_wave(2:end) +...
											 control.operators{2*k} .* sin_wave(2:end));
        end

        % build element stiffness matrix 	
		ke = zeros(dofElem,dofElem,numP,'like', control.operators{2});
		for i = 1:4
			rows = (i-1)*nRho + (1:nRho);
			for j = 1:4
				cols = (j-1)*nRho + (1:nRho);
				ke(rows,cols,:) = CE(i,j)*E + CL1(i,j)*L1 + CL2(i,j)*L2;
			end
		end

        % gradient of element stiffness matrix
        if nOutput == 3
		
			% preallocate gradient answer
			keGrad_R = cell(2*numSpin, 1);
            keGrad_L = cell(2*numSpin, 1);
		
            for k=1:numSpin

				% extend spin operators
				opX1 =  control.pwr_levels{k} * kron(CL1,control.operators{2*k-1});
				opY1 =  control.pwr_levels{k} * kron(CL1,control.operators{2*k});   				
				opX2 =  control.pwr_levels{k} * kron(CL2,control.operators{2*k-1});
				opY2 =  control.pwr_levels{k} * kron(CL2,control.operators{2*k}); 
				
				% gradient of element stiffness matrix for selected channel 
				keGrad_R{2*k-1} = kron(speye(numP), opX1);
				keGrad_R{2*k} = kron(speye(numP), opY1); 		
				keGrad_L{2*k-1} = kron(speye(numP), opX2);
				keGrad_L{2*k} = kron(speye(numP), opY2);
			
            end
			
			% combine gradient matrix
			keGrad_R = blkdiag(keGrad_R{:});      
			keGrad_L = blkdiag(keGrad_L{:});  			
        end
end

end