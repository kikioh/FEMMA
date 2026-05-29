% calculate element stiffness matrix and its gradient for solving Bloch equatoin based control problem.
% Using the 2-points linear shape function.
% Syntax:
%
%       [ke,f,keGrad] = eleMatrix_bloch(fem, control)
%
% Parameters
%
%		control.pulse_dt		- row vector, pulse time slices
%       control.offsets      	- cell array, resonance offset in Hz
%       control.pwr_level       - cell array, B1 amplitude in rad/s
%       control.waveform        - matrix, represent non-unit control waveform
% 		control.pulse_form		- string, 'xy' for Cartesian coordinate components, 'phase' for phase components
%
% Output
%
%       ke          - element stiffness matrix
%		f			- load vector in Kx=f
%       keGrad     	- gradient of element stiffness matrix wrt controls
%
% Note: the output will be used in femSolve.m to compute x and it's gradient wrt controls
%
% mengjia.he@kit.edu, 2026.04.13

function [ke,f,keGrad] = eleMatrix_bloch(fem, control)

% extract control parameters
numP = numel(control.pulse_dt);					% number of time slices
waveform = control.waveform;					% pulse waveform
dt = control.pulse_dt(1);						% time slice duration
dofNode = 3;									% dofs of each node for Bloch equatoin
nNodes = numP+1;								% number of nodes

% create bloch operators
Lx = -[0 0 0; 0 0 1; 0 -1 0];
Ly = -[0 0 -1; 0 0 0; 1 0 0];
Lz = -[0 1 0; -1 0 0; 0 0 0];
Rxyz = diag([-1/control.relax(2), -1/control.relax(2), -1/control.relax(1)]);

% coefficient for building element stiffness matrix
CE  = [-1/2  1/2; -1/2 1/2];
CL = -dt * [2/6  1/6; 1/6  2/6];

% preallocate a same-size unit matrix
E = repmat(eye(dofNode), 1, 1, numP);

% compute load vector
wf = [1; repmat(2, numP-1, 1); 1];        
f = dt/(2*control.relax(1)) * kron(wf, [0;0;1]);

% element sttiffness matrix and its gradient
switch control.pulse_form

    case 'phase'

        % transform the waveform into 3D matrix
        cos_wave = reshape(cos(waveform),1,1,[]);
        sin_wave = reshape(sin(waveform),1,1,[]);

        % addup Hamiltonian 
        L = 2*pi*control.offsets{1} * Lz + ... 
            control.pwr_levels{1} * (Lx .* cos_wave + Ly .* sin_wave) + ...
            Rxyz .* ones(1,1,numP);

        % build block-wise element stiffness matrix
        ke = zeros(2*dofNode, 2*dofNode, numP, 'like', Lx);
        for i = 1:2
            rows = (i-1)*dofNode + (1:dofNode);
            for j = 1:2
                cols = (j-1)*dofNode + (1:dofNode);
                ke(rows,cols,:) = CE(i,j)*E + CL(i,j)*L;
            end
        end

        % gradient of element stiffness matrix
        if nargout == 3

            % transform the waveform into diagonal sparse matrix
            cos_wave = control.pwr_levels{1} * spdiags(cos(waveform'), 0, numP, numP);
            sin_wave = control.pwr_levels{1} * spdiags(sin(waveform'), 0, numP, numP);

            % temporary spin oprtators
            opX =  kron(CL,Lx);
            opY =  kron(CL,Ly);

            % gradient of element stiffness matrix
            keGrad = kron(-sin_wave, opX)  + kron(cos_wave, opY);

        end

    case 'xy'

        % expansion waveform
        cos_wave = reshape(waveform(1,:), 1, 1, []);
        sin_wave = reshape(waveform(2,:), 1, 1, []);

        % addup Hamiltonian
        L = 2*pi*control.offsets{1} * Lz +...
            control.pwr_levels{1} * (Lx .* cos_wave + Ly .* sin_wave) +...
            Rxyz .*ones(1,1,numP);

        % build block-wise element stiffness matrix
        ke = zeros(2*dofNode, 2*dofNode, numP, 'like', Lx);
        for i = 1:2
            rows = (i-1)*dofNode + (1:dofNode);
            for j = 1:2
                cols = (j-1)*dofNode + (1:dofNode);
                ke(rows,cols,:) = CE(i,j)*E + CL(i,j)*L;
            end
        end

        % gradient of element stiffness matrix
        if nargout == 3

            % preallocate gradient matrix
            keGrad = cell(2, 1);

            % temporary spin oprtators
            opX =  control.pwr_levels{1} * kron(CL,Lx);
            opY =  control.pwr_levels{1} * kron(CL,Ly);

            % gradient of element stiffness matrix for selected channel
            keGrad{1} = kron(speye(numP), opX);
            keGrad{2} = kron(speye(numP), opY);

            % combine gradient matrix
            keGrad = blkdiag(keGrad{:});

        end
end

end