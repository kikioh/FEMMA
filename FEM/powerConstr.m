% calculate value and its gradient for the average power constriant in xy control. 
% Syntax:
%
%   [fval,dfdx] = powerConstr(x,Dx,pwr_constr)
%
% Parameters:
%
%   x           - numChn * N matrix, processed control variables
%
%   Dx          - 1*numChn cell array, gradient from variables precessing equation
%
%   pwr_constr  - 1*numSpin cell array, value of limited average power for each spin  
%
% Outputs:
%
%   fval        - numSpin * 1 column vector, value of power constriant for each spin
%
%   dfdx        - numSpin * nv matrix, gradient of fval to control variables
%
% mengjia.he@kit.edu, 2025.03.06


function [fval,dfdx] = powerConstr(x,pwr_constr,Dx)

% Count the outputs
nOutput = nargout;

% Check consistency
if (nOutput == 2) && (~exist('Dx','var'))
    error('gradient of variable procrssing function is required to calculate the gradient');
end

% extract control parameters
[numChn,N]  = size(x);
nv = numChn * N;
numSpin = numChn/2;

% preallocate answer
fval = zeros(numSpin,1);

% calculate constriant value and its gradient if required
if nOutput == 1

    for k = 1:numSpin

        fval(k) = (norm(x(2*k-1,:))^2 + norm(x(2*k,:))^2) / N - pwr_constr{k};

    end

elseif nOutput == 2

    dfdx = zeros(numSpin,nv);

    for k = 1:numSpin

        fval(k) = (norm(x(2*k-1,:))^2 + norm(x(2*k,:))^2) / N - pwr_constr{k};

        cols_x = 2*k-1 : numChn : (N-1)*numChn+2*k-1;
        dfdx(k,cols_x) = 1/N * 2*x(2*k-1,:) * Dx{2*k-1};

        cols_y = 2*k : numChn : (N-1)*numChn+2*k;
        dfdx(k,cols_y) = 1/N * 2*x(2*k,:) * Dx{2*k};

    end

end

end