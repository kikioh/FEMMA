% calculate the objective function, constriants and their gradient. with
% pulse shape specified in Cartesian coordinates (x and y channels). 
% Syntax:
%
%       [f0val,fval,fidelity,df0dx,dfdx] = spinSolve_xy(fem,control,xval)
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

function [f0val,fval,fidelity,df0dx,dfdx] = spinSolve_xy(fem,control,xval)

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
nv = N * numChn;
numSpin = numChn/2;

if nOutput == 3     % only function value

    % process waveform
    x = zeros(numChn,N);
    for k = 1:numChn

        % smooth variable
        xs = helmSolve(K_helm,T,xval(k,:));  

        % scale variable
        x(k,:) = sigmoid(xs);
    end

    % package control waveform
    control.waveform = x;

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
            fval = powerConstr(x,control.pwr_constr);
            
        case {'least-square','mini-1norm'}

            % calculate the objective term
            f0val =  0;

            % calculate the constraints term
            fval = zeros(n_cases*2+numSpin,1);
            fval(1:2*n_cases) = [1-fidelities; fidelities-1];
            fval(2*n_cases+1:end) = powerConstr(x,control.pwr_constr);
    end

elseif nOutput == 5                                         % function values and its gradients


    x = zeros(numChn,N); Dx = cell(1,numChn);
    for k = 1:numChn

        % smooth variable
        [xs,grad_xf] = helmSolve(K_helm,T,xval(k,:));  

        % scaled waveform and it's gradient
        [x(k,:),grad_x]  = sigmoid(xs);

        % combine smooth and scale gradient
        Dx{k} = spdiags(grad_x,0,N,N) * grad_xf;
    
    end

    % package control waveform
    control.waveform = x;
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
            [fval,dfdx] =  powerConstr(x,control.pwr_constr,Dx);   

        case {'least-square','mini-1norm'}

            % average fidelity and gradient of ensemble spin system
            [fidelities,gradients] = ensembleSpin(fem.lvn,control);

            % Add up fidelities
            fidelities=cell2mat(fidelities);
            fidelity=sum(fidelities)/n_cases;

            % objective and its gradient
            f0val =  0;
            df0dx = zeros(numChn*N,1);

            % preallocate constraint term
            fval = zeros(n_cases*2+numSpin,1);
            dfdx = zeros(n_cases*2+numSpin,nv);

            % fval and dfdx from power constriants
            [fval(2*n_cases+1:end),dfdx(2*n_cases+1:end,:)] =...
                powerConstr(x,control.pwr_constr,Dx);

            % fval from fidelity
            fval(1:2*n_cases) = [1-fidelities; fidelities-1];

            % dfdx from gradients of fidelity
            for n = 1:n_cases

                % combine LvN gradient with variable process gradient, using chain rule
                dfdx(n,:) = -combineGrad(gradients{n}, Dx);
                dfdx(n+n_cases,:) = -dfdx(n,:);

            end
    end

end