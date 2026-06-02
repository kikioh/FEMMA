% single-spin state-to-state transfer with FEMMA
%
clearvars; close all; clc;

% generate control operators
Lx = operatorCreate({'Lx'},'single','liouv');
Ly = operatorCreate({'Ly'},'single','liouv');
Lz = operatorCreate({'Lz'},'single','liouv');

% spin state
rhox = stateCreate({'Lx'},'single','liouv');
rhoz = stateCreate({'Lz'},'single','liouv');

% control parameters
control.rho_init = {rhoz};
control.rho_targ = {rhox};
control.operators = {Lx,Ly};
control.pulse_dt = 1e-6*ones(1,500);
control.pwr_levels = {2*pi*10e3*linspace(0.8,1.2,5)};
control.method = 'mma';
control.max_iter = 100;
control.tol_f = 0.995;
control.off_ops = {Lz};
control.offsets = {20e3*linspace(-0.5,0.5,51)};
control.pulse_form = 'phase';
control.integrator = 'rectangle';
control.type = 'ppt';

% create FEM system
fem.lvn = meshFEM(control,'LvN','linear');
fem.helm = meshFEM(control,'Helmholtz');

% % initial guess
guess = (2*rand(1,500)-1)*pi;
% 
% % mma parameters
[control,mma] = mmaInit(control,guess);
% 
% % run optimization
% waveform = femMMA(fem,control,mma,@spinSolve_phase);
%%
waveform = transpose(load('femma_waveform.txt'));   % 500×1 列向量，单位 rad

fidelity = fidelityCheck(control,waveform);


