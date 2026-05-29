% Calls and collect the correct amount of outputs from an objective 
% function - used by optimisation routines. Syntax:
%
%  [data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)
%
% Parameters:
%
%    x                  - objective function argument 
%    
%    objfun_handle      - handle to the objective function
%    
%    data               - data structure inherited from 
%                         fminnewton.m
%
% Outputs:
%
%    data               - modified data structure with 
%                         diagnostics from the objective
%    
%    fx                 - objective function value at x
%    
%    grad               - gradient of the objective function 
%                         at x
%    
%    hess               - Hessian of the objective function 
%                         at x
%
% Note: this function will be eliminated in a future release.
%
% david.goodwin@inano.au.dk
% i.kuprov@soton.ac.uk
%
% This function is a modified version of the original from the Spinach 
% library (MIT License). Original source: 
% <https://spindynamics.org/wiki/index.php?title=objeval.m>
%
% Modified by mengjia.he@kit.edu, 2025.11.19

function [data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)

% Check consistency
grumble(objfun_handle,x)

% Switch between function/gradient/Hessian calls
switch nargout
    
    case 2
        
        % Run the objective function for the fidelity
        [traj_data,fidelity]=feval(objfun_handle,reshape(x,data.x_shape),spin_system);
		
		if ~ismember(spin_system.control.method, {'mma', 'gcmma', 'mma-gcmma'})
        
			% Assign the data structure
			data.fx_sep_pen=fidelity;
			data.traj_data=traj_data;
			
			% Compute total error functional
			fx=fidelity(1)-sum(fidelity(2:end));
		
		else 
		
			% count number of penalty terms
			n_pen = numel(spin_system.control.p_weights);
			n_fid = numel(fidelity) - n_pen;
			fidelity_obj = transpose(cell2mat(fidelity(1:n_fid)));
			fidelity_pen = transpose(cell2mat(fidelity(n_fid+1:end)));
			
			% Assign the data structure
			data.fx_sep_pen=transpose([mean(fidelity_obj);fidelity_pen]);
			data.traj_data=traj_data;
			data.fx_square = mean((1 - fidelity_obj).^2);
			
			% Compute fx vector
			if ismember(spin_system.control.objective, {'least-square', 'mini-1norm'})
				
				fx = [1-fidelity_obj;fidelity_obj-1;fidelity_pen];
				
			elseif ismember(spin_system.control.objective, {'least-real', 'least-halfsquare'})
			
				fx = [1-fidelity_obj;fidelity_pen];
				
			else
				error('unknown objective type for MMA optimizer');
			end
		end
        
        % Increment counters
        data.count.fx=data.count.fx+1;
    
    case 3
  
		% Run the objective function for the fidelity and gradient
		[traj_data,fidelity,grad]=feval(objfun_handle,reshape(x,data.x_shape),spin_system);
			
		if ~ismember(spin_system.control.method, {'mma', 'gcmma', 'mma-gcmma'})
		
			% Assign the data structure
			data.fx_sep_pen=fidelity;
			data.traj_data=traj_data;
			
			% Compute total error functional
			fx=fidelity(1)-sum(fidelity(2:end));
			
			% Compute the total gradient
			grad=grad(:,:,1)-sum(grad(:,:,2:end),3); grad=grad(:);
			
		else
		
			% count number of penalty terms
			n_pen = numel(spin_system.control.p_weights);
			n_fid = numel(fidelity) - n_pen;
			fidelity_obj = transpose(cell2mat(fidelity(1:n_fid)));
			fidelity_pen = transpose(cell2mat(fidelity(n_fid+1:end)));
			
			% Assign the data structure
			data.fx_sep_pen=transpose([mean(fidelity_obj);fidelity_pen]);
			data.traj_data=traj_data;
			data.fx_square = mean((1 - fidelity_obj).^2);
			
			% gradient of objective and penalties
			grad_obj = transpose(cell2mat(grad(1:n_fid)));
			grad_pen = transpose(cell2mat(grad(n_fid+1:end)));
			
			% Compute fx vector and dfdx vector
			if ismember(spin_system.control.objective, {'least-square', 'mini-1norm'})
				
				fx = [1-fidelity_obj;fidelity_obj-1;fidelity_pen];
				grad = [-grad_obj; grad_obj; grad_pen];
				
			elseif ismember(spin_system.control.objective, {'least-real', 'least-halfsquare'})
			
				fx = [1-fidelity_obj;fidelity_pen];
				grad = [-grad_obj;grad_pen];
				
			else
				error('unknown objective type for MMA optimizer');
			end

		end
		
		% Increment counters
		data.count.fx=data.count.fx+1;
		data.count.gfx=data.count.gfx+1;
    
    case 4
        
        % Run the objective function for the fidelity, gradient, and Hessian
        [traj_data,fidelity,grad,hess]=feval(objfun_handle,reshape(x,data.x_shape),spin_system);
        
        % Assign the data structure
        data.fx_sep_pen=fidelity;
        data.traj_data=traj_data;
        
        % Compute total error functional
        fx=fidelity(1)-sum(fidelity(2:end));
        
        % Compute the total gradient
        grad=grad(:,:,1)-sum(grad(:,:,2:end),3); grad=grad(:);
        
        % Compute the total Hessian
        hess=hess(:,:,1)-sum(hess(:,:,2:end),3);
        
        % Increment counters
        data.count.fx=data.count.fx+1;
        data.count.gfx=data.count.gfx+1;
        data.count.hfx=data.count.hfx+1;
        
    otherwise
        
        % Complain and bomb out
        error('must have at least two output arguments.');
        
end

end

% Consistency enforcement
function grumble(objfun_handle,x)
if ~isa(objfun_handle,'function_handle')
    error('objfun_handle must be a function handle.')
end
if isempty(x)||(~isnumeric(x))||(~isreal(x))
    error('x must be a vector of real numbers')
end
end

% If you can't write it down in English, you can't code it.
%
% Peter Halpern

