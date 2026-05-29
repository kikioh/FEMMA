% Finds a local minimum of a function of several variables using Newton
% and quasi-Newton algorithms. Based on fminlbfgs.m code from D. Kroon,
% University of Twente. Syntax:
%
%          [x,data]=fminnewton(spin_system,cost_function,guess)
%
% Parameters:
%
%    spin_system   - Spinach data object that has been
%                    through optimcon.m function
%
%    cost_function - a function handle that takes the input
%                    the size of guess
%
%    guess         - the initial point of the optimisation
%
% Outputs:
%
%    x               - the final point of the optimisation
%
%    data.count.iter - iteration counter
%
%    data.count.fx   - function evaluation counter
%
%    data.count.gfx  - gradient evaluation counter
%
%    data.count.hfx  - Hessian evaluation counter
%
%    data.count.rfo  - RFO iteration counter
%
%    data.x_shape    - output of size(guess)
%
%    data.*          - further fields may be set by the
%                      objective functon
%
% david.goodwin@inano.au.dk
% i.kuprov@soton.ac.uk
%
% This function is a modified version of the original from the Spinach 
% library (MIT License). Original source:
% <https://spindynamics.org/wiki/index.php?title=fminnewton.m>
%
% Modified by mengjia.he@kit.edu, 2026.01.29

function [x,data,converg,fx_square]=fminnewton(spin_system,cost_function,guess,mma)

% extract parameters if mma is the optimization method
if ismember(spin_system.control.method, {'mma', 'gcmma', 'mma-gcmma'})
    m = mma.m;
    n = mma.n;
    xold1 = mma.xold1(:);
    xold2 = mma.xold2(:);
	scale = mma.scale;
    xmin = mma.xmin;
    xmax = mma.xmax;
    low = mma.low;
    upp = mma.upp;
    a0 = mma.a0;
    a = mma.a;
    c = mma.c;
    d = mma.d;
    kkttol = mma.kkttol;
    kktnorm = kkttol+10;

    if ~strcmp(spin_system.control.method,'mma')
        raa0 = mma.raa0;
        raa = mma.raa;
        raa0eps = mma.raa0eps;
        raaeps = mma.raaeps;
        nval = 1;
        f0valnew = 0;
		epsimin=mma.epsimin;
    end

    % Default fx and dfdx for least-square objective
    f0val = 0;
    df0dx = zeros(n,1);

    % extract the objective type
    spin_system.control.objective = mma.objective;

end

% Check consistency
grumble(spin_system,cost_function,guess);

% Initialise counters
data.count.iter=0; data.count.fx=0;
data.count.gfx=0;  data.count.hfx=0;
data.count.rfo=0;

% Look for checkpoint file
if isfield(spin_system.control,'checkpoint')&&...
        exist([spin_system.sys.scratch filesep ...
        spin_system.control.checkpoint],'file')

    % Read the initial guess from checkpoint file in global scratch
    report(spin_system,'WARNING: optimisation restarted from the checkpoint file');
    load([spin_system.sys.scratch filesep ...
        spin_system.control.checkpoint],'x');
    if numel(x)~=numel(guess)
        error('waveform size mismatch between guess and checkpoint.');
    end
    data.x_shape=size(guess);

else

    % Stretch the guess supplied by user
    data.x_shape=size(guess); x=guess(:);

end

% Prepare the freeze mask
if ~isempty(spin_system.control.freeze)

    % Validate freeze mask
    if ~all(size(spin_system.control.freeze)==data.x_shape)
        error(['freeze mask must have dimensions ' int2str(data.x_shape)]);
    end

    % Stretch freeze mask to match x
    freeze=spin_system.control.freeze(:);

else

    % Do not freeze anything
    freeze=false(size(x));

end

% Print the header
header(spin_system);

% Preallocate convergence rate
converg = zeros(1,spin_system.control.max_iter);
fx_square = zeros(1,spin_system.control.max_iter);

% Start the iteration loop
for iter=1:spin_system.control.max_iter

    % Update iteration counter
    data.count.iter=data.count.iter+1;

    % Get the search direction
    if ~ismember(spin_system.control.method, {'mma', 'gcmma', 'mma-gcmma'})
        switch spin_system.control.method

            case 'lbfgs'

                if iter==1

                    % Get objective and gradient
                    [data,fx,g]=objeval(x,cost_function,data,spin_system);

                    % Start history arrays
                    old_x=x(~freeze); dx_hist=[];
                    old_g=g(~freeze); dg_hist=[];

                    % Direction vector
                    dir=g.*(~freeze);

                    % Catch zeros
                    if (abs(fx)<1e-6)||(norm(g,2)<1e-6)
                        error('fidelity or gradient too small at iter 1, find a better guess.');
                    end

                else

                    % Update the history of dx and dg
                    dx_hist=[x(~freeze)-old_x dx_hist]; old_x=x(~freeze); %#ok<AGROW>
                    dg_hist=[g(~freeze)-old_g dg_hist]; old_g=g(~freeze); %#ok<AGROW>

                    % Truncate the history
                    if size(dx_hist,2)>spin_system.control.n_grads
                        dx_hist=dx_hist(:,1:spin_system.control.n_grads);
                    end
                    if size(dg_hist,2)>spin_system.control.n_grads
                        dg_hist=dg_hist(:,1:spin_system.control.n_grads);
                    end

                    % Get the direction
                    dir=zeros(size(freeze));
                    dir(~freeze)=-lbfgs(dx_hist,dg_hist,g(~freeze));

                end

            case {'newton','goodwin'}

                % Get objective, gradient, and Hessian
                [data,fx,g,H]=objeval(x,cost_function,data,spin_system);

                % Tidy up the inputs
                H=real(H+H')/2; g=real(g);

                % Apply freeze mask
                H=H(~freeze,~freeze);

                % Regularise the Hessian
                [H,data]=hessreg(spin_system,-H,g(~freeze),data);

                % Get the search direction
                dir=zeros(size(freeze));
                dir(~freeze)=H\g(~freeze);

        end

        % Run line search, refresh objective and gradient
        [alpha,fx,g,exitflag,data]=linesearch(spin_system,cost_function,dir,x,fx,g,data);

        % Update x
        delta_x = alpha*dir;
        x=x+delta_x;

    else

        if iter == 1

            % Get objective and gradient
			ux = x * scale;
            [data,fx,dfdx]=objeval(ux,cost_function,data,spin_system);
			dfdx = dfdx * scale;

        else

            % kkt condition
            if kktnorm <= kkttol, disp('the norm of residual vector is smaller than 0'); break; end

            % use MMA if required and current fidelity is smaller than 0.992
            if strcmp(spin_system.control.method,'mma') || (strcmp(spin_system.control.method,'mma-gcmma') && data.fx_sep_pen(1) < 0.995)

                [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,iter,x, ...
                    xmin,xmax,xold1,xold2,f0val,df0dx,fx,dfdx,low,upp,a0,a,c,d);

            else

                % use gcmma when current fidelity reaches 0.992
                % The parameters low, upp, raa0 and raa are calculated:
                [low,upp,raa0,raa] = asymp(iter,n,x,xold1,xold2,xmin,xmax,low,upp, ...
                    raa0,raa,raa0eps,raaeps,df0dx,dfdx);

                % The GCMMA subproblem is solved at the point xval:
                [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = gcmmasub(m,n,iter,epsimin,x,...
                    xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fx,dfdx,a0,a,c,d);

                % calculate function values at the point xmma
                [~,fvalnew] = objeval(xmma*scale,cost_function,data,spin_system);

                % check if the approximations are conservative
                conserv = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);

                % if approximations are non-conservative (conserv=0), repeate inner iteration
                if conserv == 0
                    for innerit = 1:16

                        % New values on the parameters raa0 and raa are calculated
                        [raa0,raa] = raaupdate(xmma,x,xmin,xmax,low,upp,f0valnew,fvalnew, ...
                            f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);

                        %%%% The GCMMA subproblem is solved with these new raa0 and raa:
                        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = gcmmasub(m,n,iter,...
                            epsimin,x,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fx,dfdx,a0,a,c,d);

                        % calculate function values at the point xmma
                        [~,fvalnew] = objeval(xmma*scale,cost_function,data,spin_system);

                        % check if the approximations have become conservative:
                        conserv = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);

                        if conserv ~= 0, break; end
                    end
                end

                % add innerit to nval
                nval = nval + innerit + 1;

            end

            % update control variables
            xold2 = xold1; xold1 = x; x = xmma;

            % Get objective and gradient
			ux = x * scale;
            [data,fx,dfdx]=objeval(ux,cost_function,data,spin_system);
			dfdx = dfdx * scale;
			
			% The residual vector of the KKT conditions is calculated:
			kktnorm = kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fx,dfdx,a0,a,c,d);

        end

        % prepare report for user
        exitflag = 0;
        alpha = 0;
        delta_x = x - xold1;
        g = transpose(mean(dfdx(1:mma.n_cases,:),1));

    end

    % Report diagnostics to user
    itrep(spin_system,g,alpha,data);

    % save convergence point
    converg(iter) = data.fx_sep_pen(1);
    if ismember(spin_system.control.method, {'mma', 'gcmma', 'mma-gcmma'})
	    fx_square(iter) = data.fx_square;
    end

    % Save checkpoint
    if isfield(spin_system.control,'checkpoint')
        save([spin_system.sys.scratch filesep spin_system.control.checkpoint],'x','-v7.3','-nocompression');
    end

    % Check termination conditions
    if (iter > 10) && (norm(delta_x,1)<spin_system.control.tol_x), exitflag = 2;

    elseif norm(g(~freeze),2)<spin_system.control.tol_g, exitflag = 1;

    elseif converg(iter)>spin_system.control.tol_f, exitflag = 3;

    end

    % Exit if necessary
    if exitflag, break; end

end

% See if iteration count was exceeded
if isempty(iter)||(iter==spin_system.control.max_iter)
    exitflag=0;
end

% Fold back the waveform
if ismember(spin_system.control.method, {'mma', 'gcmma', 'mma-gcmma'})
	x = x * scale;
end
x=reshape(x,data.x_shape);

% save convergence curve
converg = converg(1:iter);
if ismember(spin_system.control.method, {'mma', 'gcmma', 'mma-gcmma'})  
    fx_square = fx_square(1:iter);
end

% Print the footer
footer(spin_system,exitflag,data);

end

% Header printing function
function header(spin_system)
report(spin_system,'========================================================================================');
report(spin_system,'Iter  #f   #g   #H   #R    fidelity      penalties     total        alpha     |grad|    ');
report(spin_system,'----------------------------------------------------------------------------------------');
end

% Footer printing function
function footer(spin_system,exitflag,data)
report(spin_system,'----------------------------------------------------------------------------------------');
switch(spin_system.control.method)
    case 'lbfgs',  data.algorithm='LBFGS method';
    case 'newton', data.algorithm='Newton-Raphson method';
    case 'goodwin', data.algorithm='Newton-Raphson method with Goodwin acceleration';
    case {'mma','gcmma','mma-gcmma'}, data.algorithm='MMA with least-square objective';
end
switch exitflag
    case  1, message='norm(gradient,2) < tol_gfx';
    case  2, message='norm(step,1) < tol_x';
    case  0, message='number of iterations exceeded';
    case -2, message='line search found no minimum';
    case  3, message='fidelity > tol_f';
end
report(spin_system,['    Algorithm Used     : ' data.algorithm]);
report(spin_system,['    Exit message       : ' message]);
report(spin_system,['    Iterations         : ' int2str(data.count.iter)]);
report(spin_system,['    Function Count     : ' int2str(data.count.fx)]);
report(spin_system,['    Gradient Count     : ' int2str(data.count.gfx)]);
report(spin_system,['    Hessian Count      : ' int2str(data.count.hfx)]);
report(spin_system,'========================================================================================');
end

% Iteration report function
function itrep(spin_system,g,alpha,data)

% Performance figures
fid=data.fx_sep_pen(1);
pens=sum(data.fx_sep_pen(2:end));
fx = data.fx_sep_pen(1) - sum(data.fx_sep_pen(2:end));

% Apply freeze mask to gradient
if ~isempty(spin_system.control.freeze)
    g=g.*(~spin_system.control.freeze(:));
end

% Print iteration data
report(spin_system,[pad(num2str(data.count.iter,'%4.0f'),6),...
    pad(num2str(data.count.fx,'%4.0f'),5),...
    pad(num2str(data.count.gfx,'%4.0f'),5),...
    pad(num2str(data.count.hfx,'%4.0f'),5),...
    pad(num2str(data.count.rfo,'%4.0f'),5),...
    pad(num2str(fid,'%+9.6f'),11),'   '...
    pad(num2str(pens,'%+9.6f'),11),'   '...
    pad(num2str(fx,'%+9.6f'),11),'   '...
    pad(num2str(alpha,'%4.0e'),9),...
    pad(num2str(norm(g,2),'%0.3e'),10)]);

end

% Consistency enforcement
function grumble(spin_system,cost_function,guess)
if ~isa(cost_function,'function_handle')
    error('cost_function must be a function handle.');
end
if ~isfield(spin_system,'control')
    error('control field is missing from spin_system structure, run optimcon() first.');
end
if (~isnumeric(guess))||(~isreal(guess))
    error('guess must be an array of real numbers.');
end
switch spin_system.control.integrator
    case 'rectangle'
        if isempty(spin_system.control.basis)
            if size(guess,2)~=spin_system.control.pulse_nsteps
                error('the number of columns in guess must be equal to the number of time steps.');
            end
        else
            if size(spin_system.control.basis,2)~=spin_system.control.pulse_nsteps
                error('the number of columns in waveform basis must be equal to the number of time steps.');
            end
            if size(guess,2)~=size(spin_system.control.basis,1)
                error('the number of columns in guess must be equal to the number of basis functions.');
            end
        end
    case 'trapezium'
        if isempty(spin_system.control.basis)
            if size(guess,2)~=(spin_system.control.pulse_nsteps+1)
                error('the number of columns in guess must be (number of time steps)+1.');
            end
        else
            if size(spin_system.control.basis,2)~=(spin_system.control.pulse_nsteps+1)
                error('the number of columns in waveform basis must be (number of time steps)+1.');
            end
            if size(guess,2)~=size(spin_system.control.basis,1)
                error('the number of columns in guess must be equal to the number of basis functions.');
            end
        end
    otherwise
        error('unknown time propagation algorithm.');
end
end

% There is a sacred horror about everything grand. It is easy to
% admire mediocrity and hills; but whatever is too lofty, a geni-
% us as well as a mountain, an assembly as well as a masterpiece,
% seen too near, is appalling... People have a strange feeling of
% aversion to anything grand. They see abysses, they do not see
% sublimity; they see the monster, they do not see the prodigy.
%
% Victor Hugo - Ninety-three

