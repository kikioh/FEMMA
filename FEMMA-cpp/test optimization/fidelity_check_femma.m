% Fidelity check of a FEMMA-exported pulse against the MATLAB reference.
%
% Parses the exported Bruker shape file, reads the offset/power/duration grids
% from its header into the control struct, then recomputes the fidelity.  Only
% the spin-system fields (operators, states, type, method) are set by hand; the
% sampling grids come straight from the pulse file, so the two never drift.

clearvars; close all; clc;

% generate control operators
Lx = operatorCreate({'Lx'},'single','liouv');
Ly = operatorCreate({'Ly'},'single','liouv');
Lz = operatorCreate({'Lz'},'single','liouv');

% spin state
rhox = stateCreate({'Lx'},'single','liouv');
rhoz = stateCreate({'Lz'},'single','liouv');

% spin-system / method fields (NOT stored in the shape file)
control.rho_init   = {rhoz};
control.rho_targ   = {rhox};
control.operators  = {Lx,Ly};
control.off_ops    = {Lz};
control.method     = 'mma';
control.max_iter   = 100;
control.tol_f      = 0.995;
control.integrator = 'rectangle';
control.type       = 'ppt';

% ---- parse the exported shape file ---------------------------------------
shapefile = 'D:\Cpp\femma-cpp\femma_pulse.txt';
[amp,pha,hdr] = parseShape(shapefile);

% pull offsets, powers, pulse_dt, and pulse_form from the header
control = controlFromShape(control,hdr);

% report what was read back from the pulse
fprintf('read from %s:\n',shapefile);
fprintf('  points   : %d\n',          numel(pha));
fprintf('  duration : %s us\n',        hdr.FEMMA_DURATION_US);
fprintf('  offsets  : %s kHz\n',       hdr.FEMMA_OFFSET_RANGE_KHZ);
fprintf('  power    : %s%% of %s kHz\n',hdr.FEMMA_PWR_RANGE_PCT,hdr.FEMMA_PWR_NOMINAL_KHZ);
fprintf('  form     : %s\n',           control.pulse_form);

% ---- waveform and fidelity -----------------------------------------------
% phase form: amplitude constant, only the phase matters; degrees -> radians
waveform = (pha.') * pi/180;

fidelity = fidelityCheck(control,waveform);
fprintf('\nfidelity stored in file = %s\n',hdr.FEMMA_FIDELITY);
fprintf('fidelity recomputed     = %.6f\n',fidelity);
