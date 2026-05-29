% check the fidelity of optimal pulse 
%
% mengjia.he@kit.edu, 2025.03.12


function fidelity = fidelityCheck(control,waveform,spin_system)

% Default spin system
if ~exist('spin_system','var')
    % 11.74T magnet
    sys.magnet = 11.74;
    sys.isotopes = {'1H'};
    inter.zeeman.scalar={0};
    
    % Basis set
    bas.formalism='sphten-liouv';
    bas.approximation='none';
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

end

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

% pulse slices
pulse_dt = control.pulse_dt;

% setup initial state
rhoz = state(spin_system,'Lz','1H');
rhoz = rhoz/norm(rhoz,'fro');
dim = size(control.rho_init{1},1);

% transform phase into xy controls
switch control.pulse_form
    case 'phase'
        func = @(k) {cos(waveform(k,:)), sin(waveform(k,:))}; 
        numSpin = size(waveform,1);
    case 'xy'
        func = @(k) {waveform(2*k-1,:), waveform(2*k,:)};     
        numSpin = size(waveform,1)/2;
end

% package waveform into cell array
wave = cell(1, 2*numSpin);  
for k = 1:numSpin, wave(2*k-1:2*k) = func(k);  end

% Preallocate fidelity
fidelity = zeros(1,n_cases);

% Loop over total cases
parfor n=1:n_cases  %#ok<*PFBNS>

    % Extract ensemble indices
    n_rho=catalog(n,1); n_pwr=catalog(n,2); n_off=catalog(n,3);

    % Get initial and target state
    rho_init=control.rho_init{n_rho};
    rho_targ=control.rho_targ{n_rho};

    % Preallocate local pulse and Hamiltonian
    H = sparse(size(control.operators{1},1),size(control.operators{1},2));
    local_pulse = cell(size(wave));

    % Loop over spins
    for k = 1:numSpin

        % local pulse
        local_pulse{2*k-1} = control.pwr_levels{k}(n_pwr)* wave{2*k-1};
        local_pulse{2*k} = control.pwr_levels{k}(n_pwr)* wave{2*k};

        % local Hamiltonian
        H = H + 2*pi*control.offsets{k}(offset_index{n_off}(k))*control.off_ops{k};
    end

    % calculate fidelity
    switch  control.type
        case 'ppt'

			% Spinach v2.8
			fidelity(n) =  rho_targ'*shaped_pulse_xy(spin_system,H,control.operators,local_pulse,pulse_dt,rho_init);
			
			% Spinach v2.10 and later
			% fidelity(n) =  rho_targ'*shaped_pulse_xy(spin_system,H,control.operators,local_pulse,pulse_dt,rho_init,'expv-pwc');

        case 'gate'
		
			% Spinach v2.8
			[~,~,Prop] =  shaped_pulse_xy(spin_system,H,control.operators,local_pulse,pulse_dt,rhoz);
			
			% Spinach v2.10 and later
            % [~,~,Prop] =  shaped_pulse_xy(spin_system,H,control.operators,local_pulse,pulse_dt,rhoz,'expv-pwc');
			
			
            fidelity(n) =real(trace(rho_targ'*Prop))/dim;
    end
end

% average fidelity
fidelity = mean(fidelity);

end