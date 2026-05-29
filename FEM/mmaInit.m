% Validate optimal control options and prepare the MMA parameters. Syntax:
%
%       [mma,control] = mmaInit(control,guess)
%
% Parameters:
%
%       control         - control data structure, specify the optimal control parameters
%
%       guess           - initial waveform
%
% Outputs:
%
%       control         - updated control structure
%
%       mma             - mma data structure
%
% Note: this program is modified from https://spindynamics.org/wiki/index.php?title=Optimcon.m
%
% mengjia.he@kit.edu, 2025.02.04

function [control,mma] = mmaInit(control,guess)

% Default termination conditions
if ~isfield(control, 'tol_f'), control.tol_f = 0.999; end

% Default objective type
if ~isfield(control, 'objective'), control.objective = 'least-square'; end

% Default maximum iterations
if ~isfield(control, 'max_iter'), control.max_iter = 100; end

% Default pulse form
if ~isfield(control, 'pulse_form'), control.pulse_form = 'xy'; end

% Process control type
if ~isfield(control,'type'), control.type='ppt'; end

% Process control type
if ~isfield(control,'equation'), control.equation='LvN'; end


% Default control method
if ~isfield(control, 'method')
    control.method = 'mma';
    disp('The default optimizatoin method mma is used.');
end


% Special process for LvN equation
if strcmp(control.equation,'LvN')

	% Process control operators and power levels
	if ~isfield(control,'operators') || ~isfield(control,'pwr_levels')

		% Complain and bomb out
		error('control.operators and control.pwr_levels must be supplied.');

	else

		% Input validation
		if (~iscell(control.operators))||(~all(cellfun(@ismatrix,control.operators(:))))
			error('control.operators must be a cell array of matrices.');
		end

		% Basic type validation
		if ~iscell(control.pwr_levels)
			error('control.pwr_levels must be a cell array of row vectors.');
		end
		if numel(control.pwr_levels)~=numel(control.operators)/2
			error('control.pwr_levels and control.pwr_levels must be aligned.');
		end

		% Content validation
		for n=1:numel(control.pwr_levels)
			if (~isnumeric(control.pwr_levels{n}))||...
					(~isreal(control.pwr_levels{n}))||...
					(~isrow(control.pwr_levels{n}))
				error('elements of control.pwr_levels must be row vectors of real numbers.');
			end
		end

		pwr_level_lens = cellfun(@length, control.pwr_levels);
		if numel(unique(pwr_level_lens)) ~= 1
			error('the power level for each spin should be aligned.')
		end
		

	end
	
	% Process drift generators
	if isfield(control, 'drifts')
		
		if ~isnumeric(control.drifts)
			error('control.drifts must be a numeric matrix');
		end
		
		if (size(control.drifts,1) ~= size(control.operators{1},1)) || (size(control.drifts,2) ~= size(control.operators{1},2))
			error('drift Hamiltonian must have the same dimension with control.operators');
		end
			
		if ismatrix(control.drifts)
			control.drifts = control.drifts .* ones(1,1,size(guess,2));
		end
		
		if ndims(control.drifts) ==  3
		
			if size(control.drifts,3) ~= size(guess,2)
				error('time-dependent drift Hamiltonians must be aligned with guess');
			end
			
		end
		
	else
		
		control.drifts = zeros(size(control.operators{2},1), size(control.operators{2},2), size(guess,2), 'like', control.operators{2});
		fprintf('drift Hamiltonian not supplied, use all-zero drift Hamiltonian\n');
	end

	% Process initial waveform 
	if strcmp(control.pulse_form,'phase') && size(guess,1) ~= numel(control.operators)/2
		error('the number of columns in guess must be equal to the number of control operators.');
	end
	
	% Process offset distributions
	if isfield(control,'offsets')&&isfield(control,'off_ops')

		% Basic type validation
		if ~iscell(control.offsets)
			error('control.offsets must be a cell array of row vectors.');
		end
		if ~iscell(control.off_ops)
			error('control.offsets must be a cell array of square matrices.');
		end
		if numel(control.offsets)~=numel(control.off_ops)
			error('control.offsets and control.off_ops must have the same size.');
		end

		% Content validation
		for n=1:numel(control.offsets)
			if (~isnumeric(control.offsets{n}))||...
					(~isreal(control.offsets{n}))||...
					(~isrow(control.offsets{n}))
				error('elements of control.offsets must be row vectors of real numbers.');
			end
			if (~isnumeric(control.off_ops{n}))||...
					(size(control.off_ops{n},1)~=size(control.off_ops{n},2))
				error('elements of control.off_ops must be square matrices.');
			end
		end

	else

		% Default is no offsets
		control.offsets={};
		control.off_ops={};

	end

end

% number of power levels 
n_power_levls = numel(control.pwr_levels{1});

% number of offsets
n_offsets=prod(cellfun(@numel,control.offsets));

% Process power constriant
if strcmp(control.pulse_form,'phase') && isfield(control, 'pwr_constr')

    disp('power constriants will be ignored in amplitude-phase control');

elseif strcmp(control.pulse_form,'xy')

    if ~isfield(control, 'pwr_constr')

        % Default power constriant values
        control.pwr_constr = num2cell(2*ones(1,numel(control.pwr_levels)));

    else
        % Content validation
        for n=1:numel(control.pwr_constr)
            if (~isnumeric(control.pwr_constr{n}))||...
                    (~isreal(control.pwr_constr{n}))||...
                    ((control.pwr_constr{n})<=0)
                error('elements of control.pwr_constr must be positive numbers, a good start is 2.');
            end
        end
    end

end

% Process initial states
if isfield(control,'rho_init')

    % Input validation
    if ~iscell(control.rho_init)
        error('control.rho_init must be a cell array of vectors.');
    end
    if strcmp(control.type,'ppt')
        for n=1:numel(control.rho_init)
            if (~isnumeric(control.rho_init{n}))||(~iscolumn(control.rho_init{n}))
                error('control.rho_init must be a cell array of column vectors.');
            end
        end
    end
else

    % Complain and bomb out
    error('initial states must be supplied in control.rho_init field.');

end

% Process target states
if isfield(control,'rho_targ')

    % Input validation
    if ~iscell(control.rho_targ)
        error('control.rho_targ must be a cell array of vectors.');
    end
    if strcmp(control.type,'ppt')
        for n=1:numel(control.rho_targ)
            if (~isnumeric(control.rho_targ{n}))||(~iscolumn(control.rho_targ{n}))
                error('control.rho_targ must be a cell array of column vectors.');
            end
        end
    end
    if numel(control.rho_targ)~=numel(control.rho_init)
        error('control.rho_targ must have the same size as control.rho_init');
    end
    n_state_pairs=numel(control.rho_init);
else

    % Complain and bomb out
    error('target states must be supplied in control.rho_targ field.');

end

% Process pulse form
if ~isfield(control,'pulse_form'), error('the pulse form should be specified with ap or xy'); end

% Process initial waveform
if (~isnumeric(guess))||(~isreal(guess))
    error('guess must be an array of real numbers.');
end
if size(guess,2) ~= numel(control.pulse_dt)
    error('the column size of initial guess should be equal to number of time slices.');
end
if strcmp(control.pulse_form,'xy') && size(guess,1) ~= numel(control.operators)
    error('the number of columns in guess must be equal to the number of control operators.');
end

[numChn,numP] = size(guess);
control.numChn = numChn;

% ensemble size
control.n_cases = n_offsets * n_power_levls * n_state_pairs;
mma.n_cases = control.n_cases;

% Setup MMA parameters
switch control.pulse_form

    case 'xy'
        numSpin = numChn/2;
        mma.n = numChn*numP;
        mma.xmin = -ones(mma.n,1);
        mma.xmax = ones(mma.n,1);
		mma.scale = 1;

    case 'phase'
        numSpin = numChn;
        mma.n = numChn*numP;
		
		% bound used for different optimization in optimal control paper
		switch control.type

			case {'ppt','ppt-auxmat','ppt-FD'} % point-to-point transfer
				mma.scale = 1;
				mma.xmin = -15*pi * ones(mma.n,1);
				mma.xmax = 15*pi * ones(mma.n,1);

			% because the the variable processing (helmholtz filter) in the optimization; if scale is futher applied, the 
			% optimization may become highly non-linear, so when combined with variable processing in FEM, the suggest scale is 1.
			% in this setting, MMA can get stuck in a local optimum even for a target fidelity of 0.995.
			
			case 'gate' % propagator optimization
				mma.scale = 10*pi;
				mma.xmin = -ones(mma.n,1);
				mma.xmax = ones(mma.n,1);	
			% here auxmat with exact propagation is used to compute gradient and Hessian, so no variable processing is needed, variable 
			% scale helps MMA to converg more stablely.
		end
		
		
end
		
mma.epsimin = 1e-7;
mma.xval = guess;
mma.xold1 = mma.xval;
mma.xold2 = mma.xval;
mma.low = mma.xmin;
mma.upp = mma.xmax;
mma.maxoutit  = control.max_iter;
mma.kkttol  = 0;
control.scale = mma.scale;

% setup a, b, c, d parameters
mma.objective = control.objective;
switch control.objective

    % case    'mini-real'          % real part of overlap

        % mma.m = numSpin;
        % mma.a0 = 1;
        % mma.a = zeros(numSpin,1);
        % mma.c = 1000 * ones(numSpin,1);
        % mma.d = ones(numSpin,1);

        % if ~strcmp(control.method,'mma')
            % mma.raa0    = 1e-2;
            % mma.raa     = 1e-2 * ones(numSpin,1);
            % mma.raa0eps = 1e-6;
            % mma.raaeps  = 1e-6 * ones(numSpin,1);
        % end

    case 'least-square'           % least-square objective of the overlap, each fidelity cooresponds to 2 constraints
        p = control.n_cases;
        mma.m = 2*p+numSpin;
        mma.a0 = 1;
        mma.a = zeros(mma.m,1);
        mma.c = [zeros(2*p,1);1000 * ones(numSpin,1)];
        mma.d = [2*ones(2*p,1);zeros(numSpin,1)];

        if ~strcmp(control.method,'mma')
            mma.raa0    = 1e-2;
            mma.raa     = 1e-2 * ones(mma.m,1);
            mma.raa0eps = 1e-5;
            mma.raaeps  = 1e-5 * ones(mma.m,1);
        end

    case 'mini-1norm'           % minimum 1–norm of the overlap, each fidelity cooresponds to 2 constraints
        p = control.n_cases;
        mma.m = 2*p+numSpin;
        mma.a0 = 1;
        mma.a = zeros(mma.m,1);
        mma.c = [ones(2*p,1);1000 * ones(numSpin,1)];
        mma.d = zeros(mma.m,1);

        if ~strcmp(control.method,'mma')
            mma.raa0    = 1e-2;
            mma.raa     = 1e-2 * ones(mma.m,1);
            mma.raa0eps = 1e-6;
            mma.raaeps  = 1e-6 * ones(mma.m,1);
        end
	
    case 'least-halfsquare'     % least-square objective of the overlap, each fidelity cooresponds to a constraint
        p = control.n_cases;
        mma.m = p+numSpin;
        mma.a0 = 0.01;
        mma.a = zeros(mma.m,1);
        mma.c = [zeros(p,1);1000 * ones(numSpin,1)];
        mma.d = [2*ones(p,1);zeros(numSpin,1)];

        if ~strcmp(control.method,'mma')
            mma.raa0    = 1e-2;
            mma.raa     = 1e-2 * ones(mma.m,1);
            mma.raa0eps = 1e-5;
            mma.raaeps  = 1e-5 * ones(mma.m,1);
        end

    case 'least-real'     		% least-square objective of the overlap, each fidelity cooresponds to a constraint
        p = control.n_cases;
        mma.m = p+numSpin;
        mma.a0 = 1;
        mma.a = zeros(mma.m,1);
        mma.c = [ones(p,1); 1000 * ones(numSpin,1)];
        mma.d = zeros(mma.m,1);

        if ~strcmp(control.method,'mma')
            mma.raa0    = 1e-2;
            mma.raa     = 1e-2 * ones(mma.m,1);
            mma.raa0eps = 1e-5;
            mma.raaeps  = 1e-5 * ones(mma.m,1);
        end
end

end