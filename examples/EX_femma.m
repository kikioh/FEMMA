% benchmark convergence rate and speed of FEM for 1-spin optimal
% contorl
%
% mengjia.he@kit.edu, 2025.02.18

clearvars; close all; clc;

% % generate control operators
Lx = operatorCreate({'Lx'},'single','liouv');
Ly = operatorCreate({'Ly'},'single','liouv');
Lz = operatorCreate({'Lz'},'single','liouv');

% specify spin state
rhox = stateCreate({'Lx'},'single','liouv');
rhoz = stateCreate({'Lz'},'single','liouv');

% create control
control.pulse_dt = 1e-6*ones(1,500);                    % time duration
control.rho_init = {rhoz};                              % initial state
control.rho_targ = {rhox};                              % target state
control.operators = {Lx,Ly};                            % control operators
control.off_ops = {Lz};                                 % offset operators
control.offsets = {20e3*linspace(-0.5,0.5,51)};         % Offset distribution, Hz
control.pwr_levels = {2*pi*10e3*linspace(0.8,1.2,5)};   % Pulse power ensemble, rad/s
control.max_iter = 50;                                  % Termination condition
control.tol_f = 0.995;                                  % Minimum fidelity
control.method = 'mma';                                 % optimization method
control.pulse_form = 'phase';                           % control pulse form
control.integrator = 'rectangle';                       % approximation type
control.type = 'ppt';

% create fem system
fem.lvn = meshFEM(control,'LvN','linear');              % FEM parameters for LvN equation
fem.helm = meshFEM(control,'Helmholtz');                % FEM parameters for Helmholtz

% initial guess
guess = pi*(2*rand(1,numel(control.pulse_dt))-1);

% create mma parameters
[control,mma] = mmaInit(control,guess);

% run optimization
[waveform,converg] = femMMA(fem,control,mma,@spinSolve_phase);

% check fidelity
% fid = fidelityCheck(control,waveform);
% fprintf('Absolute fidelity difference between FEM and Spinach: %.6f\n',...
% abs(fid - converg(end)));