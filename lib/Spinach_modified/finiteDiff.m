% Gradient Ascent Pulse Engineering (GRAPE) objective function and gradient.
% Propagates the system through a user-supplied shaped pulse
% from a given initial state and projects the result onto the given final
% state. The fidelity is returned, along with its gradient and Hessian 
% with respect to amplitudes of all control operators at every time step
% of the shaped pulse. 
%
% Using second-order center finite difference to compute the gradient of step 
% propogator wrt the controls, the step size was chosen to be h=1e-1.
% ref: 2011_JMR_Second order gradient ascent pulse engineering, 10.1016/j.jmr.2011.07.023
% 
% Syntax:
%
%  [traj_data,fidelity,grad]=finiteDiff(spin_system,drifts,controls,...
%                                       waveform,rho_init,rho_targ,...
%                                       fidelity_type)
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
% david.goodwin@inano.au.dk
% uluk.rasulov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% This function is a modified version of the original from the Spinach 
% library (MIT License). Original source:
% <https://spindynamics.org/wiki/index.php?title=grape.m>
%
% Modified by mengjia.he@kit.edu, 2026.01.17

function [fwd_traj,fidelity,grad]=finiteDiff(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type) %#ok<*PFBNS>
    
% Count the outputs
n_outputs=nargout();

% Extract the timing grid
dt=spin_system.control.pulse_dt;

% Run array preallocations
switch spin_system.control.integrator

    % Piecewise-constant
    case 'rectangle'

        % Number of time intervals and control operators
        nsteps=size(waveform,2); nctrls=size(waveform,1);

        % Preallocate forward and backward trajectories
        fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);
        bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);

		% Parfor needs these empty declarations
		fwd_dP={}; bwd_dP={}; fwd_d2P={}; P={};

    % Piecewise-linear
    case 'trapezium'
      
        % Number of time intervals and control operators
        nsteps=size(waveform,2)-1; nctrls=size(waveform,1);

        % Preallocate forward and backward trajectories
        fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);
        bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i);

        % Cannot do trapezium Hessians yet
        fwd_dP={}; bwd_dP={}; fwd_d2P={}; P={};
        
    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Hush up the output
spin_system.sys.output='hush';

% Initialise forward and backward trajectories
fwd_traj(:,1)=rho_init; bwd_traj(:,1)=rho_targ;

% Parallelisation prep
if n_outputs>2

    % Strip the spin system object for communication efficiency
    ss_parfor.sys=spin_system.sys; ss_parfor.tols=spin_system.tols;
    ss_parfor.bas.formalism=spin_system.bas.formalism;

    % Define a vector and a matrix of zeroes for auxiliary systems
    zero_state=spalloc(size(rho_init,1),size(rho_init,2),0);
    zero_drift=spalloc(size(drifts{1},1),size(drifts{1},2),0);

end

% Run forward and backward propagation
switch spin_system.control.integrator

    % Piecewise-constant
    case 'rectangle'

        % Loop over time steps
        for n=1:nsteps

            % Decide current drifts
            if numel(drifts)==1

                % Time-independent drifts, including
                % conj-transpose for dissipative cases
                L_forw=drifts{1}; L_back=drifts{1}';

            else

                % Time-dependent drifts, including
                % conj-transpose for dissipative cases
                L_forw=drifts{n}; L_back=drifts{nsteps+1-n}';

            end

            % Add current controls to current drifts, including
            % conjugate-transpose for dissipative controls; the
            % waveform is always real
            for k=1:nctrls
                L_forw=L_forw+waveform(k,n)*controls{k};
                L_back=L_back+waveform(k,nsteps+1-n)*controls{k}';
            end

            % Take a time step forwards and backwards
            fwd_traj(:,n+1)=step(spin_system,L_forw,fwd_traj(:,n),+dt(n));
            bwd_traj(:,n+1)=step(spin_system,L_back,bwd_traj(:,n),-dt(nsteps+1-n));

        end

    % Piecewise-linear
    case 'trapezium'

        % Loop over time steps
        for n=1:nsteps

            % Decide current drifts
            if numel(drifts)==1

                % Time-independent drifts, including
                % conjugate-transpose for dissipative cases
                L_forw_left=drifts{1};  L_forw_right=drifts{1};
                L_back_left=drifts{1}'; L_back_right=drifts{1}';

            else

                % Time-dependent drifts, including
                % conjugate-transpose for dissipative cases
                L_forw_left=drifts{n};           L_forw_right=drifts{n+1};
                L_back_left=drifts{nsteps+1-n}'; L_back_right=drifts{nsteps+2-n}';

            end

            % Add current controls to current drifts, including
            % conjugate-transpose for dissipative controls; the
            % waveform is always real
            for k=1:nctrls
                L_forw_left=L_forw_left+waveform(k,n)*controls{k};
                L_forw_right=L_forw_right+waveform(k,n+1)*controls{k};
                L_back_left=L_back_left+waveform(k,nsteps+1-n)*controls{k}';
                L_back_right=L_back_right+waveform(k,nsteps+2-n)*controls{k}';
            end

            % Take a time step forwards and backwards
            fwd_traj(:,n+1)=step(spin_system,{ L_forw_left,...
                                              (L_forw_left+L_forw_right)/2,...
                                               L_forw_right},fwd_traj(:,n),+dt(n));
            bwd_traj(:,n+1)=step(spin_system,{ L_back_right,...
                                              (L_back_right+L_back_left)/2,...
                                               L_back_left},bwd_traj(:,n),-dt(nsteps+1-n));
            
        end

    otherwise

        % Complain and bomb out
        error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

end

% Calculate the state overlap
overlap=rho_targ'*fwd_traj(:,end);

% Compute gradient
if n_outputs>2

	% Set finite difference step
	h = 1e-1;
	
    % Flip the backward trajectory
    bwd_traj=fliplr(bwd_traj);
    
    % Preallocate results
    grad=zeros(size(waveform),'like',1i);
    
    % Integrator-specific paths
    switch spin_system.control.integrator

        % Piecewise-constant
        case 'rectangle'

            % Loop over control sequence
            for n=1:nsteps
                
                % Allocate local gradient column
                grad_col=zeros(nctrls,1,'like',1i);

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
				
				% Pick the forward trajectory state
				aux_vec=fwd_traj(:,n);
				
                % Calculate gradient at this timestep
                for k=1:nctrls
            
                    % Add the tiny step to Hamiltonian		 
					L_fwd = L + h*controls{k};
					L_bwd = L - h*controls{k};

                    % Propagate the auxiliary vector
                    aux_vec_fwd=step(ss_parfor,L_fwd,aux_vec,dt(n));
					aux_vec_bwd=step(ss_parfor,L_bwd,aux_vec,dt(n));
            
                    % Compute the derivative
                    grad_col(k)=bwd_traj(:,n+1)'*(aux_vec_fwd - aux_vec_bwd) / (2*h);
            
                end
        
                % Add to the gradient array
                grad(:,n)=grad_col;

            end

        % Piecewise-linear
        case 'trapezium'

            % Loop over control sequence
            for n=1:(nsteps+1)

                % Allocate local gradient column
                grad_col=zeros(nctrls,1,'like',1i);

                % First step is special
                if n==1

                    % Decide current drift
                    if numel(drifts)==1

                        % Time-independent drift
                        L={drifts{1},drifts{1}};

                    else

                        % Time-dependent drift
                        L={drifts{1},drifts{2}};

                    end
    
                    % Loop over controls
                    for k=1:nctrls

                        % Build the auxiliary matrix
                        [G_fwd_l,G_bwd_l,~,~]=generator_trap(L,controls,dt(n),waveform(:,1),waveform(:,2),k,h);

                        % Build the auxiliary vector
                        aux_vec=fwd_traj(:,n);

                        % Propagate the auxiliary vector
                        aux_vec_fwd=step(ss_parfor,G_fwd_l,aux_vec,dt(n));
						aux_vec_bwd=step(ss_parfor,G_bwd_l,aux_vec,dt(n));

                        % Compute the derivative
                        grad_col(k)=bwd_traj(:,n+1)'*(aux_vec_fwd - aux_vec_bwd) / (2*h);

                    end
                
                % Last step is special
                elseif n==(nsteps+1)

                    % Decide current drift
                    if numel(drifts)==1

                        % Time-independent drift
                        L={drifts{1},drifts{1}};

                    else

                        % Time-dependent drift
                        L={drifts{n-1},drifts{n}};

                    end

                    % Loop over controls
                    for k=1:nctrls
                        
                        % Build the auxiliary matrix
                        [~,~,G_fwd_r,G_bwd_r]=generator_trap(L,controls,dt(n-1),waveform(:,(end-1)),waveform(:,end),k,h);
                        
                        % Build the auxiliary vector
                        aux_vec=fwd_traj(:,n-1);

                        % Propagate the auxiliary vector
                        aux_vec_fwd=step(ss_parfor,G_fwd_r,aux_vec,dt(n-1));
						aux_vec_bwd=step(ss_parfor,G_bwd_r,aux_vec,dt(n-1));
						
                        % Compute the derivative
                        grad_col(k)=bwd_traj(:,end)'*(aux_vec_fwd - aux_vec_bwd) / (2*h);

                    end

                % Middle steps
                else

                    % Loop over controls
                    for k=1:nctrls

                        % Left pair of drifts
                        if numel(drifts)==1

                            % Time-independent drift
                            L={drifts{1},drifts{1}};

                        else

                            % Time-dependent drift
                            L={drifts{n},drifts{n+1}};

                        end
                        
                        % Auxiliary matrix
                        [G_fwd_l,G_bwd_l,~,~]=generator_trap(L,controls,dt(n),waveform(:,n),waveform(:,n+1),k,h);

                        % Build the auxiliary vector
                        aux_vec_a= fwd_traj(:,n);

						% Propagate the auxiliary vector
                        aux_vec_fwd=step(ss_parfor,G_fwd_l,aux_vec_a,dt(n));
						aux_vec_bwd=step(ss_parfor,G_bwd_l,aux_vec_a,dt(n));

                        % Product rule: [dP2]*[P1]*rho part
                        grad_col(k)=grad_col(k)+bwd_traj(:,n+1)'*(aux_vec_fwd - aux_vec_bwd) / (2*h);
                        
                        % Right pair of drifts
                        if numel(drifts)==1

                            % Time-independent drift
                            L={drifts{1},drifts{1}};

                        else

                            % Time-dependent drift
                            L={drifts{n-1},drifts{n}};

                        end

                        % Auxiliary vector and matrix
                        [~,~,G_fwd_r,G_bwd_r]=generator_trap(L,controls,dt(n-1),waveform(:,n-1),waveform(:,n),k,h);

                        % Build the auxiliary vector
                        aux_vec_b= fwd_traj(:,n-1);
                        
                        % Propagate the auxiliary vector
                        aux_vec_fwd=step(ss_parfor,G_fwd_r,aux_vec_b,dt(n-1));
						aux_vec_bwd=step(ss_parfor,G_bwd_r,aux_vec_b,dt(n-1));

                        % Product rule: [P2]*[dP1]*rho part
                        grad_col(k)=grad_col(k)+bwd_traj(:,n)'*(aux_vec_fwd - aux_vec_bwd) / (2*h);

                    end

                end

                % Add to the gradient array
                grad(:,n)=grad_col;

            end

        otherwise

            % Complain and bomb out
            error('unknown integrator: must be ''rectangle'' or ''trapezium''.');

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

end
