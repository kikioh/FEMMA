% Main optimization function that calling MMA algorithm to solve the
% quantum spin optimal control. Syntax:
%
%   [waveform,fid,nval] = femMMA(fem,control,mma,objHandle)
%
% Parameters:
%
%   fem         - primary FEM data structure, created by meshFEM.m
%
%   control     - control data structure, specify the optimal control
%               parameters
%
%   mma         - mma data structure
%
%   objHandle   - spinSolve_xy or spinSolve_phase, calculate the objective
%               terms and constriant terms
%
% Outputs:
%
%   waveform    - optimized pulse waveform
%
%   fid         - convergence of fidelity
%
%   nval        - number of inner iterations for GCMMA method
%
% mengjia.he@kit.edu, 2025.02.04

function [waveform,fid,xval,nval] = femMMA(fem,control,mma,objHandle)

% extract mma parameters
m = mma.m;
n = mma.n;
epsimin = mma.epsimin;
xval = mma.xval(:);
xold1 = mma.xold1(:);
xold2 = mma.xold2(:);
xmin = mma.xmin;
xmax = mma.xmax;
low = mma.low;
upp = mma.upp;
a0 = mma.a0;
a = mma.a;
c = mma.c;
d = mma.d;
kkttol = mma.kkttol;
maxoutit = mma.maxoutit;
kktnorm = kkttol+10;
if ~strcmp(control.method,'mma')
    raa0 = mma.raa0;
    raa = mma.raa;
    raa0eps = mma.raa0eps;
    raaeps = mma.raaeps;
end
numChn = control.numChn;

% calculate function values and gradients for initial guess
[f0val,fval,fidelity,df0dx,dfdx] = objHandle(fem,control,mma.xval);


% Preallocate convergence answer
fid = zeros(1,control.max_iter);   % save the fidelity as outeriter increase
outeriter = 1;                     % outer iteration number
nval = 1;                          % times of solving system equations
innerit = 0;                       % Default inner iteration

% inform the user
fprintf('=============================================\n');
fprintf('Iter      Inner     Fidelity       |grad|    \n');
fprintf('---------------------------------------------\n');

% current gradient amplitude
switch control.objective
    case 'mini-real', grad = norm(df0dx,2);
    case {'least-square','mini-1norm'}, grad = mean(vecnorm(dfdx(1:control.n_cases,:),2,2));
end

% report current fidelity and gradient
fprintf([pad(num2str(outeriter,'%4.0f'),10),...
         pad(num2str(innerit,'%4.0f'),10),...
         pad(num2str(fidelity,'%0.5f'),15),...
         pad(num2str(grad,'%0.3e'),10) '\n']);

% Check termination conditions
if  fidelity>=control.tol_f
    disp('target fidelity achieved');
    fid = fidelity;
else

    % record fidelity of initial guess
    fid(1) = fidelity;

    % start optimization
    for outeriter = 1:maxoutit

        % break loop if kktnorm is smaller than 0
        if kktnorm <= kkttol, disp('the norm of residual vector is smaller than 0'); break; end

        % use MMA if required and current fidelity is smaller than 0.992
        if strcmp(control.method,'mma') || (strcmp(control.method,'mma-gcmma') && fidelity < 0.995)

            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,outeriter,xval, ...
                xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);

        else

            % use gcmma when current fidelity reaches 0.992
            % The parameters low, upp, raa0 and raa are calculated:
            [low,upp,raa0,raa] = asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
                raa0,raa,raa0eps,raaeps,df0dx,dfdx);

            % The GCMMA subproblem is solved at the point xval:
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = gcmmasub(m,n,outeriter,epsimin,xval,...
                xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);

            % calculate function values at the point xmma
            [f0valnew,fvalnew,~] = objHandle(fem,control,reshape(xmma,numChn,[]));

            % check if the approximations are conservative
            conserv = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);

            % if approximations are non-conservative (conserv=0), repeate inner iteration
            if conserv == 0
                for innerit = 1:16

                    % New values on the parameters raa0 and raa are calculated
                    [raa0,raa] = raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
                        f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);

                    %%%% The GCMMA subproblem is solved with these new raa0 and raa:
                    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = gcmmasub(m,n,outeriter,...
                        epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);

                    % calculate function values at the point xmma
                    [f0valnew,fvalnew,~] = objHandle(fem,control,reshape(xmma,numChn,[]));

                    % check if the approximations have become conservative:
                    conserv = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);

                    if conserv ~= 0, break; end
                end
            end

            % add innerit to nval
            nval = nval + innerit + 1;

        end

        %%%% No more inner iterations. Some vectors are updated:
        xold2 = xold1;
        xold1 = xval;
        xval  = xmma;

        % calculate function values and gradients at xval.
        [f0val,fval,fidelity,df0dx,dfdx] = objHandle(fem,control,reshape(xval,numChn,[]));

        % The residual vector of the KKT conditions is calculated:
        kktnorm = kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);

        % record fidelity
        fid(outeriter+1) = fidelity;

        % current gradient amplitude
        switch control.objective
            case 'mini-real', grad = norm(df0dx,2);
            case {'least-square','mini-1norm'}, grad = mean(vecnorm(dfdx(1:control.n_cases,:),2,2));
        end

        % report current fidelity and gradient
        fprintf([pad(num2str(outeriter+1,'%4.0f'),10),...
            pad(num2str(innerit,'%4.0f'),10),...
            pad(num2str(fidelity,'%0.5f'),15),...
            pad(num2str(grad,'%0.3e'),10) '\n']);

        % Check termination conditions
        if  fidelity>=control.tol_f, disp('target fidelity achieved'); break; end

    end

    % save convergence curve
    fid = fid(1:outeriter+1);
end

% process waveform
T = sum(control.pulse_dt);
xval = reshape(xval,numChn,[]);
waveform = zeros(numChn,numel(control.pulse_dt));
for k = 1:numChn
    if strcmp(control.pulse_form, 'xy')
        xs = helmSolve(fem.helm.K, T, xval(k,:));
        waveform(k,:) = sigmoid(xs);
    elseif strcmp(control.pulse_form, 'phase')
        waveform(k,:) = helmSolve(fem.helm.K, T, xval(k,:), 1, 0);
		waveform(k,:) = waveform(k,:) * control.scale;
        % waveform(k,:) = mod(waveform(k,:), 2*pi);
    end
end

end

