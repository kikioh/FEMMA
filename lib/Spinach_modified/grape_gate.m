% Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
% and Hessian. Propagates the system through a user-supplied shaped pulse
% from a given initial state and projects the result onto the given final
% state. The fidelity is returned, along with its gradient and Hessian
% with respect to amplitudes of all control operators at every time step
% of the shaped pulse. This function is for the Hilbert space version of GRAPE 
% for propagator optimization. Syntax:
%
% [traj_data,fidelity,grad]=grape_gate(spin_system,drifts,controls,waveform,init,targ,fidelity_type)
%
% Parameters:
%
%   spin_system         - Spinach data object that has been through
%                         the optimcon.m problem setup function.
%
%   drifts              - the drift Liouvillians: a cell array con-
%                         taining one matrix (for time-independent
%                         drift) or multiple matrices (for time-de-
%                         pendent drift).
%
%   controls            - control operators in Liouville space (cell
%                         array of matrices).
%
%   waveform            - control coefficients for each control ope-
%                         rator (in rows of a matrix), rad/s
%
%   rho_init            - initial state of the system as a vector in
%                         Liouville space.
%
%   rho_targ            - target state of the system as a vector in
%                         Liouville space.
%
%   fidelity_type       - 'real'   (real part of the overlap)
%                         'imag'   (imaginary part of the overlap)
%                         'square' (absolute square of the overlap)
%
% Outputs:
%
%   fidelity            - fidelity of the control sequence
%
%   grad                - gradient of the fidelity with respect to
%                         the control sequence
%
%   traj_data.forward   - forward trajectory from the initial condi-
%                         tion(a stack of state vectors)
%
%   traj_data.backward  - backward trajectory from the target state
%                         (a stack of state vectors)
%
% Note: this is a low level function that is not designed to be called
%       directly. Use grape_xy.m and grape_phase.m instead.
%
% david.goodwin@inano.au.dk
% uluk.rasulov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% This function is a modified version of the original from the Spinach 
% library (MIT License). Original source:
% <https://spindynamics.org/wiki/index.php?title=grape.m> 
%
% Modified by mengjia.he@kit.edu, 2026.01.27

function [traj_data,fidelity,grad]=grape_gate(spin_system,drifts,controls,waveform,init,targ,fidelity_type) %#ok<*PFBNS>

% Check consistency
grumble(spin_system,drifts,controls,waveform,init,targ);

% Count the outputs
n_outputs=nargout();

% Extract the timing grid
dt=spin_system.control.pulse_dt;

% Number of time intervals and control operators
nsteps=size(waveform,2); nctrls=size(waveform,1);

% matrix dimension of hamiltonian
dim=size(drifts{1},1); dim_sz = [dim,dim]; nrm = dim;

% Preallocate forward and backward propagators
fwd_traj = init; bwd_traj = targ;

% Hush up the output
spin_system.sys.output='hush';

% Parallelisation prep
if n_outputs>2

    % Strip the spin system object for communication efficiency
    ss_parfor.sys=spin_system.sys; ss_parfor.tols=spin_system.tols;
    ss_parfor.bas.formalism=spin_system.bas.formalism;

    % Define a vector and a matrix of zeroes for auxiliary systems
    zero_drift=spalloc(size(drifts{1},1),size(drifts{1},2),0);

    % Preallocate results
    grad=zeros(size(waveform),'like',1i);

    % initialise propagator cell
    dP_rho=cell(1,nctrls); dP_rho(:)={zeros([dim_sz,nsteps])};

end

% Run forward and backward propagation

if n_outputs<=2

    % Loop over time steps
    for n=1:nsteps

        % Decide current drifts
        if numel(drifts)==1

            % Time-independent drifts, including
            % conj-transpose for dissipative cases
            L=drifts{1};

        else

            % Time-dependent drifts, including
            % conj-transpose for dissipative cases
            L=drifts{n};

        end

        % Add current controls to current drifts, including
        % conjugate-transpose for dissipative controls; the
        % waveform is always real
        for k=1:nctrls
            L=L+waveform(k,n)*controls{k};
        end

        % Take a time step forwards and backwards
        fwd_traj=step(spin_system,L,fwd_traj,+dt(n));

    end

    % Calculate the state overlap
    overlap=trace(targ'*fwd_traj)/nrm;


else % Compute gradient

    % Loop over control sequence
    for n=1:nsteps

        % Decide the current drift
        if numel(drifts)==1

            % Time-independent drift
            L=drifts{1};

        else

            % Time-dependent drift
            L=drifts{n};

        end

        % Add current controls
        for k=1:nctrls
            L=L+waveform(k,n)*controls{k};
        end

        % Create auxiliary state vector
        aux_vec=[zero_drift; fwd_traj];

        % Calculate gradient at this timestep
        for k=1:nctrls

            % Create auxiliary system
            aux_matrix=[ L           controls{k}
                zero_drift  L           ];

            % Propagate the auxiliary vector
            aux_vector=step(ss_parfor,aux_matrix,aux_vec,dt(n));

            % Compute the derivative
            dP_rho{k}(:,:,n)=aux_vector(1:(end/2),:);

        end

        % Take a time step forwards and backwards
        fwd_traj=step(spin_system,L,fwd_traj,+dt(n));

    end

    % Calculate the state overlap
    overlap=trace(targ'*fwd_traj)/nrm;

    % loop backwards over timesteps (for efficiency)
    for n = nsteps:-1:1

        % Decide the current drift
        if numel(drifts)==1

            % Time-independent drift
            L=drifts{1};

        else

            % Time-dependent drift
            L=drifts{n};

        end

        % calculate the gradient elements
        for k=1:nctrls

            grad(k,n) = trace(bwd_traj'*dP_rho{k}(:,:,n))/nrm;

            % add the current controls
            L=L+waveform(k,n)*controls{k};
        end

        % adjoint state
        bwd_traj=step(spin_system,L,bwd_traj,-dt(n));
    end

end

% Fidelity and its derivatives
switch fidelity_type

    case {'real'}

        % Real part of the overlap
        fidelity=real(overlap);

        % Update gradient
        if exist('grad','var'), grad=real(grad); end

    case {'imag'}

        % Imaginary part of the overlap
        fidelity=imag(overlap);

        % Update gradient
        if exist('grad','var'), grad=imag(grad); end

    case {'square'}

        % Absolute square of the overalp
        fidelity=overlap*conj(overlap);

        % Update gradient
        if exist('grad','var')

            % Product rule
            grad=grad*conj(overlap)+overlap*conj(grad);

            % Cleaning up
            grad=real(grad);

        end

    otherwise

        % Complain and bomb out
        error('unknown fidelity type');

end

% Return trajectory data
traj_data=[];

% Catch unreachable objectives
if abs(fidelity)==0
    spin_system.sys.output=1;
    report(spin_system,'exactly zero fidelity: either the target is unreachable');
    report(spin_system,'from the source, or the initial guess is very poor.');
    error('GRAPE cannot proceed.');
end
if exist('grad','var')&&(norm(grad,1)==0)
    spin_system.sys.output=1;
    report(spin_system,'exactly zero gradient: either the target is unreachable');
    report(spin_system,'from the source, or the initial guess is very poor.');
    error('GRAPE cannot proceed.');
end

end

% Consistency enforcement
function grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('optimal control module requires Lioville space formalism.');
end
if (~isnumeric(rho_init))
    error('rho_init must be a matrix.');
end
if (~isnumeric(rho_targ))
    error('rho_targ must be a matrix.');
end
if ~iscell(drifts)
    error('drifts must be a cell array of matrices.');
end
for n=1:numel(drifts)
    if (~isnumeric(drifts{n}))||(size(drifts{n},1)~=size(drifts{n},2))
        error('all elements of drifts cell array must be square matrices.');
    end
    if (size(drifts{n},1)~=size(rho_init,1))||...
            (size(drifts{n},1)~=size(rho_targ,1))
        error('dimensions of drift, rho_init and rho_targ must be consistent.');
    end
end
if ~iscell(controls)
    error('controls must be a cell array of square matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||...
            (size(controls{n},1)~=size(controls{n},2))||...
            (size(controls{n},1)~=size(drifts{1},1))
        error('control operators must have the same size as drift operators.');
    end
end
if (~isnumeric(waveform))||(~isreal(waveform))
    error('waveform must be a real numeric array.');
end
if size(waveform,1)~=numel(controls)
    error('number of waveform rows must be equal to the number of controls.');
end
if size(waveform,2)~=spin_system.control.pulse_ntpts
    error(['waveform must have ' int2str(spin_system.control.pulse_ntpts) ' columns.']);
end
if strcmp(spin_system.control.integrator,'rectangle')&&...
        (size(spin_system.control.pulse_dt,2)~=size(waveform,2))
    error('pulse_dt must have the same length as waveform for rectangle integrator');
end
if strcmp(spin_system.control.integrator,'trapezium')&&...
        (size(spin_system.control.pulse_dt,2)+1~=size(waveform,2))
    error('pulse_dt must be one element shorter than waveform for trapezium integrator');
end
end

% In any culture, subculture, or family in which belief is valued above
% thought, self-surrender is valued above self-expression, and conformity
% is valued above integrity, those who preserve their self-esteem are
% likely to be heroic exceptions.
%
% Nathaniel Branden

