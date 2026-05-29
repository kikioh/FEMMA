% A parallel wrapper that enables ensemble optimal control optimisations. 
% This function handles systems with multiple control power levels, multiple 
% resonance offsets. Syntax:
%
%      [fidelity,gradient]=ensembleSpin(fem,control)
%
% Parameters:
%
% 	fem		- FEM data structure
%
%	control	- control data structure
%
% Output
%
%   fidelity   - averaged fidelity for ensemble spin system
%
%   gradient   - gradient of final spin state to control waveform
%
% Note: this funciton has refered to a Spinach ensemble.m function:
% <https://spindynamics.org/wiki/index.php?title=ensemble.m>
%
% mengjia.he@kit.edu, 2025.02.18

function [fidelities, gradients] = ensembleSpin(fem,control)

% number od Output
numOut = nargout;

% offset ensemble size
off_ens_sizes=cellfun(@numel,control.offsets);
numSpin = numel(off_ens_sizes);
if ~isempty(off_ens_sizes)

    % Multi-channel offsets index
    n_offset_vals=prod(off_ens_sizes);
    rev_lens = off_ens_sizes(end:-1:1);                     % Reverse the order of the length array
    rev_stride = [1, cumprod(rev_lens(1:end-1))];           % Compute reversed strides (row-major order)
    stride = rev_stride(end:-1:1);                          % Reverse strides back to the original order
    offset_index = cell(1,n_offset_vals);
    for k = 1:n_offset_vals

        % Column-major order index calculation (the last array changes the fastest)
        offset_index{k} = mod(floor((k-1) ./ stride), off_ens_sizes) + 1;
    end
else
    n_offset_vals=1;
end

% Extract ensemble grid dimensions
n_state_pairs=numel(control.rho_init);     	            % State-target pair count
n_power_levls=numel(control.pwr_levels{1});             % Power level count

% Create a catalog of the ensemble
catalog=(1:n_state_pairs)';
catalog=[kron(ones(n_power_levls,1),catalog) kron((1:n_power_levls)',ones(size(catalog,1),1))];
catalog=[kron(ones(n_offset_vals,1),catalog) kron((1:n_offset_vals)',ones(size(catalog,1),1))];

% count total number of cases
n_cases=size(catalog,1);

% extract initial state and waveform
offsets = control.offsets;
pwr_levels = control.pwr_levels;
dofNode = numel(control.rho_init{1});

% Preallocate outputs
fidelities=cell(n_cases,1); gradients=cell(n_cases,1);

% % number of workers
nworkers=poolsize;

% for n=1:n_cases
parfor (n=1:n_cases,nworkers) %#ok<*PFBNS>

    % Extract ensemble indices
    n_rho=catalog(n,1); n_pwr=catalog(n,2); n_off=catalog(n,3);

    % extract control parameters
    ctmp  = control;
	
    % Get initial and target state
	ctmp.rho_init=control.rho_init(n_rho);
    ctmp.rho_targ=control.rho_targ(n_rho);
	rho_targ = control.rho_targ{n_rho};
	
	% resonance offsets and power levels
    ctmp.offsets = cell(1,numSpin);
    for k = 1:numSpin
        ctmp.pwr_levels{k} = pwr_levels{k}(n_pwr);
        ctmp.offsets{k} = offsets{k}(offset_index{n_off}(k));
    end

    if numOut == 1

        % stiffness of LvN equation
        [K,f] = globalStiff(fem,ctmp);

        % solution of the linear equations
        rho_traj = femSolve(ctmp,K,f);

        % fidelity
		switch control.integrator
		
			case 'rectangle', fidelities{n} = real(rho_targ' * rho_traj(end-dofNode+1:end));
				
			case 'trapezium', fidelities{n} = real(rho_targ' * rho_traj(end-dofNode*2+1:end-dofNode));
			
		end

    elseif numOut == 2

		switch control.integrator
		
			case 'rectangle'
			
				% stiffness of LvN equation
				[K, f, Kgrad] = globalStiff(fem,ctmp);
			
				% solution of the linear equations
				[rho_traj, gradients{n}] = femSolve(ctmp,K,f,Kgrad);
				
				% fidelity
				fidelities{n} = real(rho_targ' * rho_traj(end-dofNode+1:end));
		
		case 'trapezium' 
		
				% stiffness of LvN equation
				[K, f, Kgrad_L, Kgrad_R] = globalStiff(fem, ctmp);
			
				% solution of the linear equations
				[rho_traj, gradients{n}] = femSolve(ctmp, K, f, Kgrad_L, Kgrad_R);
				
				% fidelity
				fidelities{n} = real(rho_targ' * rho_traj(end-dofNode*2+1:end-dofNode));
				
		end
		
    end

end

end

