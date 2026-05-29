% calculate the objective function, constriants and their gradient, the
% pulse shapes are given by phases and the amplitudes are fixed at 1.
% Syntax:
%
%   [f0val,fval,fidelity,df0dx,dfdx] = spinSolve_phase(fem,control,xval)
%
% Parameters:
%
% 	fem		    - FEM data structure
%
%	control	    - control data structure
%
%   xval        - the current control variables
%
% Output
%
%   f0val       - objective for MMA optimization
%
%   df0dx       - gradient of f0val to xval
%
%   fidelity    - average fidelity over total cases
%
%   fval        - constriants for MMA optimization
%
%   dfdx        - gradient of fval to xval
%
% mengjia.he@kit.edu, 2025.02.18

function [f0val,fval,fidelity,df0dx,dfdx] = spinSolve_phase(fem,control,xval)

% Count the outputs
nOutput = nargout;

% count the number of ensemble spin system
n_cases = control.n_cases;             

% extract fem parameters
K_helm = fem.helm.K;

% extract control parameters
T = sum(control.pulse_dt);
N = numel(control.pulse_dt);
numChn = control.numChn;
numSpin = numChn;
nv = numChn * N;
scale = control.scale;
% objWeight = repmat(conj(control.rho_targ{1}),1,nv);

if nOutput == 3                         % only function value
    
    % smoothed waveform
    xs = zeros(numChn,N);
    for k = 1:numChn
        xs(k,:) = helmSolve(K_helm,T,xval(k,:),1,0);  
		
    end
	
    % update control waveform
    control.waveform = scale * xs;

    % average fidelity of ensemble spin system
    fidelities = ensembleSpin(fem.lvn,control);

    % Add up fidelities
    fidelities=cell2mat(fidelities);
    fidelity=sum(fidelities)/n_cases;

    switch control.objective

        case 'mini-real'
            % calculate the objective term
            f0val =  1-fidelity;

            % calculate the constraints term
            fval = zeros(numSpin,1);

        case {'least-square','mini-1norm'}

            % calculate the objective term
            f0val =  0;

            % calculate the constraints term
            fval = zeros(n_cases*2+numSpin,1);
            fval(1:2*n_cases) = [1-fidelities; fidelities-1];
            fval(n_cases*2+1:end) = 0;
    end

elseif nOutput == 5                    % function values and its gradients

	% smoothed waveform
    xs = zeros(numChn,N); Dx = cell(1,numChn);
    for k = 1:numChn

        % smoothed pulse shape
        [xs(k,:),Dx{k}] = helmSolve(K_helm,T,xval(k,:),1,0);  
		Dx{k} = Dx{k} * scale;
    
    end

    % update control waveform
    control.waveform = scale * xs;

    switch control.objective

        case 'mini-real'

            % average fidelity and gradient of ensemble spin system
            [fidelities,gradients] = ensembleSpin(fem.lvn,control);

            % Add up fidelities
            fidelities=cell2mat(fidelities);
            fidelity=sum(fidelities)/n_cases;

            % Add up gradients
            drho_dxs = zeros(size(gradients{1}));
            for n = 1:n_cases
                drho_dxs = drho_dxs + gradients{n};
            end
            drho_dxs = drho_dxs/n_cases;
            
            % calculate the objective term
            f0val =  1-fidelity;

            % combine LvN gradient with variable process gradient, using chain rule
            df0dx = -transpose(combineGrad(drho_dxs, Dx));

            % calculate the constraints term
            fval = zeros(numSpin,1);
            dfdx = zeros(numSpin,nv);

        case {'least-square','mini-1norm'}

            % average fidelity and gradient of ensemble spin system
            [fidelities,gradients] = ensembleSpin(fem.lvn,control);

            % Add up fidelities
            fidelities=cell2mat(fidelities);
            fidelity=sum(fidelities)/n_cases;

            % objective and its gradient
            f0val =  0;
            df0dx = zeros(nv,1);

            % value of constraints
            fval = zeros(n_cases*2+numSpin,1);
            fval(1:2*n_cases) = [1-fidelities; fidelities-1];
            fval(n_cases*2+1:end) = 0;

            % prelocate gradient of constraints
            dfdx = zeros(n_cases*2+numSpin,nv);

            % gradient of energy constraints
            dfdx(n_cases*2+1:end,:) = 0;

            % Add up gradients
            for n = 1:n_cases

                % combine LvN gradient with variable process gradient, using chain rule
                dfdx(n,:) = -combineGrad(gradients{n}, Dx);
                dfdx(n+n_cases,:) = -dfdx(n,:);

            end

    end

end