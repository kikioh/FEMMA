% benchmark convergence rate of lbfgs/newton/mma (specified in control.method) 
% for single-spin propagator optimization
%
% use GRAPE and auxiliary matrix formalism to compute the gradient and Hessian
%
% mengjia.he@kit.edu, 2025.02.21

% clear variables
clearvars; close all; clc;

% specify experiment parameters
control.pulse_dt = 1e-6*ones(1,500);                  % time duration
control.offsets = {20e3*linspace(-0.5,0.5,51)};       % Offset distribution, Hz
control.pwr_levels = {2*pi*10e3*linspace(0.9,1.1,5)}; % Pulse power ensemble, rad/s

% spin system
sys.magnet = 1;
sys.isotopes = {'1H'};
inter.zeeman.scalar={0};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the control operators
Lx=operator(spin_system,'Lx','1H');
Ly=operator(spin_system,'Ly','1H');
Lz=operator(spin_system,'Lz','1H');

% Calculate the Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% setup target propagator
P_init = speye(size(Lx));
P_targ = propagator(spin_system,Lx,pi/2);
control.rho_init={P_init};                                      % initial propagator
control.rho_targ={P_targ};                                      % target propagator

% Define control parameters
control.type = 'gate';											% propagator optimization
control.pulse_form = 'phase';									% waveform type
control.drifts={{H}};                                           % Drifts
control.off_ops = {Lz};                                         % resonance offset opeators
control.operators = {Lx,Ly};                                    % Controls
control.method='lbfgs';                                         % Optimisation method
control.max_iter=100;                                           % Maximum iterations
control.tol_f=0.995;                                          	% Maximum fidelity
% control.penalties={'DNS'};                                      % Penalty types
% control.p_weights=1e-9;                                         % Penalty weights
control.parallel='ensemble';                                    % Parallelisation mode
control.freeze=zeros(1,numel(control.pulse_dt));                % Freeze mask
control.amplitudes=ones(1,numel(control.pulse_dt));             % Amplitude profile

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% initial guess
guess= pi/3 * (2*rand(1,numel(control.pulse_dt))-1);

% run optimization with lbfgs/newton
pulse_shape=fminnewton(spin_system,@grape_phase,guess);

% % run optimization with mma
% [control,mma] = mmaInit(control,guess);
% pulse_shape=fminnewton(spin_system,@grape_phase,guess,mma);