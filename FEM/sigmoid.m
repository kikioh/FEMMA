% scale function for limiting the control variables to [-1, 1]
%
% mengjia.he@kit.edu, 2025.01.22

function [x,grad_x] = sigmoid(xf,k2)

% Count the outputs
n_outputs=nargout;

% Default exp coefficient
if ~exist('k2','var'), k2 = 10; end

% scaled waveform and it's gradient
if n_outputs == 1
	x = -1+2./(1+exp(-k2*xf));
	
elseif n_outputs == 2 
	
	x = -1+2./(1+exp(-k2*xf));
    xf = transpose(xf);
	grad_x = 2*k2*exp(-k2*xf)./(1+exp(-k2*xf)).^2;
    
	
end
end