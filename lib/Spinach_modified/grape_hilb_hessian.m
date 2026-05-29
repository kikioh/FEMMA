% Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
% and Hessian. Propagates the system through a user-supplied shaped pulse
% from a given initial state and projects the result onto the given final
% state. The fidelity is returned, along with its gradient and Hessian
% with respect to amplitudes of all control operators at every time step
% of the shaped pulse. This function is for the Hilbert space version of GRAPE.
%
% Syntax:
%
%  [traj_data,fidelity,grad,hess]=grape_hilb_hessian(spin_system,drifts,controls,waveform,init,targ,fidelity_type)
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
% library (MIT License). Original source: <https://spindynamics.org/wiki/index.php?title=grape.m>
%
% Modified by mengjia.he@kit.edu, 2025.11.07

function [traj_data,fidelity,grad,hess]=grape_hilb_hessian(spin_system,drifts,controls,waveform,init,targ,fidelity_type) %#ok<*PFBNS>

% Check consistency
grumble(spin_system,drifts,controls,waveform,init,targ);

% compute the effective propagator in Liouville space, as 
% U_eff = U_N * U_N-1 * ... * U_2 * U_1 applys to both Liouville space and Hilbert space
if strcmp(spin_system.bas.formalism, 'zeeman-hilb')

spin_system.bas.formalism = 'zeeman-liouv';

end

% Count the outputs
n_outputs=nargout();

% Extract the timing grid
dt=spin_system.control.pulse_dt;

% Number of time intervals and control operators
nsteps=size(waveform,2); nctrls=size(waveform,1);

% matrix dimension of hamiltonian
dim=size(drifts{1},1); dim_sz = [dim,dim]; nrm = dim;

if n_outputs<=3

	% Preallocate forward and backward propagators
	fwd_traj = init; bwd_traj = targ;

else 

	% Preallocate forward and backward trajectories
	fwd_traj=zeros([size(init,1) size(init,2) (nsteps+1)],'like',1i);
	bwd_traj=zeros([size(init,1) size(init,2) (nsteps+1)],'like',1i);

	% Preallocate arrays used in Hessian calculation
	% Both Newton and Goodwin
	fwd_dP=cell(nctrls,nsteps);
	bwd_dP=cell(nctrls,nsteps); 
	fwd_d2P=cell(nctrls,nctrls,nsteps);

	% Goodwin method precomputes propagators
	P=cell(1,nsteps);
	
	% Initialise forward and backward trajectories
	fwd_traj(:,:,1)=init; bwd_traj(:,:,1)=targ;

end

% Hush up the output
spin_system.sys.output='hush';

% Parallelisation prep
if n_outputs==3

    % Strip the spin system object for communication efficiency
    ss_parfor.sys=spin_system.sys; ss_parfor.tols=spin_system.tols;
    ss_parfor.bas.formalism=spin_system.bas.formalism;

    % Define a vector and a matrix of zeroes for auxiliary systems
    zero_drift=spalloc(size(drifts{1},1),size(drifts{1},2),0);

    % Preallocate results
    grad=zeros(size(waveform),'like',1i);

    % initialise propagator cell
    dP_rho=cell(1,nctrls); dP_rho(:)={zeros([dim_sz,nsteps])};
	
elseif n_outputs==4

	% Strip the spin system object for communication efficiency
    ss_parfor.sys=spin_system.sys; ss_parfor.tols=spin_system.tols;
    ss_parfor.bas.formalism=spin_system.bas.formalism;

    % Define a vector and a matrix of zeroes for auxiliary systems
    zero_drift=spalloc(size(drifts{1},1),size(drifts{1},2),0);

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


elseif n_outputs == 3 % Compute gradient

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

else % cumpute Hessian 

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

		% Goodwin's Hessian method pre-computes cumulative propagators
		P{n}=propagator(ss_parfor,L_forw,dt(n));
		if n>1
			P{n}=P{n}*P{n-1};
			P{n}=clean_up(spin_system,P{n},spin_system.tols.prop_chop);
		end

		% Take a time step forwards and backwards
		fwd_traj(:,:,n+1)=step(spin_system,L_forw,fwd_traj(:,:,n),+dt(n));
		bwd_traj(:,:,n+1)=step(spin_system,L_back,bwd_traj(:,:,n),-dt(nsteps+1-n));

	end

	% Calculate the state overlap
	overlap=trace(targ'*fwd_traj(:,:,end))/nrm;
	
	% Compute gradient
	% Flip the backward trajectory
    bwd_traj=flip(bwd_traj,3);
    
    % Preallocate results
    grad=zeros(size(waveform),'like',1i);
	
	% Loop over control sequence
	parfor n=1:nsteps
		
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

		% Calculate gradient at this timestep
		for k=1:nctrls
	
			% Create auxiliary system
			aux_matrix=[ L           controls{k}
						 zero_drift  L           ];

			% Build the auxiliary vector
			aux_vec=[zero_drift; fwd_traj(:,:,n)];
	
			% Propagate the auxiliary vector
			aux_vec=step(ss_parfor,aux_matrix,aux_vec,dt(n));
	
			% Compute the derivative
			grad_col(k)=trace(bwd_traj(:,:,n+1)'*aux_vec(1:(end/2),:))/nrm;
	
		end

		% Add to the gradient array
		grad(:,n)=grad_col;

	end
	
	% compute Hessian
	 % Flip the backward trajectory for ease of paralellisation
    bwd_traj=flip(bwd_traj,3);

    % Compute derivative trajectories
    parfor n=1:nsteps

        % Decide the current drift
        if numel(drifts)==1

            % Time-independent drifts, including conj-transpose for dissipative cases
            L_forw=drifts{1}; L_back=drifts{1}';

        else

            % Time-dependent drifts, including conj-transpose for dissipative cases
            L_forw=drifts{n}; L_back=drifts{nsteps+1-n}';

        end

        % Add current controls to current drifts, including conjugate-transpose 
		% for dissipative controls; the waveform is always real
        for k=1:nctrls
            L_forw=L_forw+waveform(k,n)*controls{k};
            L_back=L_back+waveform(k,nsteps+1-n)*controls{k}';
        end

        % Preallocate local arrays
        fwd_dP_col=cell(nctrls,1);
        bwd_dP_col=cell(nctrls,1);
        fwd_d2P_block=cell(nctrls,nctrls);

        % Loop over control pairs
        for k=1:nctrls
            for j=1:nctrls

                % Diagonal elements
                if k==j

                    % Create forward auxiliary matrix
                    aux_matrix=[L_forw        controls{k}   zero_drift
                                zero_drift    L_forw        controls{j}
                                zero_drift    zero_drift    L_forw];

                    % Create forward auxiliary vector
                    aux_vec=[zero_drift; zero_drift; fwd_traj(:,:,n)];

                    % Propagate the auxiliary vector
                    aux_vec=step(spin_system,aux_matrix,aux_vec,dt(n));

                    % Only store acton by dP on rho
                    fwd_dP_col{k}=aux_vec((end/3+1):(2*end/3),:);

                    % Only store action by d2P on rho
                    fwd_d2P_block{k,j}=2*aux_vec(1:(end/3),:);

                    % Create backward auxiliary matrix
                    aux_matrix=[L_back      controls{k}'
                                zero_drift  L_back     ];

                    % Create backward auxiliary vector
                    aux_vec=[zero_drift; bwd_traj(:,:,n)];

                    % Propagate the auxiliary vector
                    aux_vec=step(spin_system,aux_matrix,aux_vec,-dt(nsteps+1-n));

                    % Only store the action by dP on rho
                    bwd_dP_col{k}=aux_vec(1:(end/2),:);

                % Off-diagonal elements
                else

                    % Create the forward auxiliary matrix
                    aux_matrix=[L_forw        controls{k}   zero_drift
                                zero_drift    L_forw        controls{j}
                                zero_drift    zero_drift    L_forw];

                    % Create forward auxiliary vector
                    aux_vec=[zero_drift; zero_drift; fwd_traj(:,:,n)];

                    % Propagate the auxiliary vector
                    aux_vec=step(spin_system,aux_matrix,aux_vec,dt(n));

                    % Only store action by d2P on rho
                    fwd_d2P_block{k,j}=2*aux_vec(1:(end/3),:);

                end

            end
        end

        % Store derivative trajectories
        fwd_dP(:,n)=fwd_dP_col; 
        bwd_dP(:,n)=bwd_dP_col;
        fwd_d2P(:,:,n)=fwd_d2P_block;

    end

    % Preallocate Hessian matrix
    hess=zeros(nctrls,nsteps,nctrls,nsteps,'like',1i);

    % Off-diagonal Hessian elements
	% Flip the backwards derivative trajectory
	bwd_dP=fliplr(bwd_dP);

	% Loop over timesteps
	parfor n=1:nsteps

		% Loop over controls
		for k=1:nctrls

			% Propagate forward derivatives to first time step
			fwd_dP{k,n}=P{n}'*fwd_dP{k,n};

			% From second step
			if n>1

				% Propagate backward derivatives to first time step
				bwd_dP{k,n}=bwd_dP{k,n}'*P{n-1};

			end

		end

	end

	% Flip the backward trajectory
	bwd_traj=flip(bwd_traj,3);

	% Loop over timesteps
	parfor n=1:nsteps

		% Preallocate local Hessian column
		hess_col=zeros(nctrls,nsteps,nctrls,1,'like',1i);

		% Outer control loop
		for k=1:numel(controls)

			% From second step
			if n>1

				% Inner control loop
				for j=1:numel(controls)

					% Construct array of forward derivatives
					% array_fwd_dP=cat(2,fwd_dP{j,1:n-1});
					array_fwd_dP=cat(3,fwd_dP{j,1:n-1});
					
					% Multiply out current backward derivatives and
					% array of all forward derivatives
					% hess_col(j,1:n-1,k)=bwd_dP{k,n}*array_fwd_dP;	
					hess_col(j,1:n-1,k)=transpose(squeeze(sum(bwd_dP{k,n}.*permute(array_fwd_dP,[2 1 3]),[1 2])))/nrm;

				end

			end

			% Inner control loop
			for j=1:numel(controls)

				% Calculate non-mixed derivatives
				hess_col(j,n,k)=trace(bwd_traj(:,:,n+1)'*fwd_d2P{k,j,n})/nrm;

			end

		end

		% Add to Hessian array
		hess(:,:,:,n)=hess_col;

	end

	% Merge the blocks and reorder
    hess=reshape(hess,nsteps*nctrls,nsteps*nctrls);
    
    % Force Hessian symmetry and fill empty Hessian entries
    hess=(hess+hess.').*(~kron(eye(nsteps),ones(nctrls,nctrls)))+...
         (hess+hess.').*( kron(eye(nsteps),ones(nctrls,nctrls)))/2;

end

% Fidelity and its derivatives
switch fidelity_type

    case {'real'}

        % Real part of the overlap
        fidelity=real(overlap);

        % Update gradient
        if exist('grad','var'), grad=real(grad); end
		
		% Update Hessian
        if exist('hess','var'), hess=real(hess); end

    case {'imag'}

        % Imaginary part of the overlap
        fidelity=imag(overlap);

        % Update gradient
        if exist('grad','var'), grad=imag(grad); end
		
		% Update Hessian
        if exist('hess','var'), hess=imag(hess); end

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
		
		% Update Hessian
        if exist('hess','var')
            
            % Product rule
            hess=hess*conj(overlap)+grad(:)*transpose(conj(grad(:)))+...
                 conj(grad(:))*transpose(grad(:))+conj(hess)*overlap;
            
            % Cleaning up
            hess=real(hess);
        
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
% if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    % error('optimal control module requires Lioville space formalism.');
% end
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

