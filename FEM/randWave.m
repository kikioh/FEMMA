% creat random pulse shape for FEM test. Syntax:
%
% 		[wave,wave_xy_pl] = randWave(fem.helm, control)
%
% Parameters:
%
%	fem.helm	- fem structure of helmholtz filter
%	control		- optimal control structure
%
% Output
%
%	wave 	- a matrix represents the pulse shape , the size depends on the pulse form and integrator type
%
%	wave_xy_pl 	- wave in Cartesian coordinates with power level weighted, a matrix represents the pulse shape,
%				  the size depends on the pulse form and integrator type
%
% mengjia.he@kit.edu, 2025.10.16

function [wave, wave_xy, wave_xy_pl] = randWave(helm, control)

% extract parameters
K = helm.K;							% stiffness matrix for Helmholtz filter
T = sum(control.pulse_dt);			% time duration of pulse
numP = numel(control.pulse_dt); 	% number of time slices
numChn = numel(control.operators); 	% number of contorl channels

% generate random pulse shape
switch control.pulse_form

    case 'phase'
        switch control.integrator
            case 'rectangle'
                wave = zeros(numChn/2, numP);
                for m =1:numChn/2
                    % xc = pi*(2*rand(1,numP)-1);
					
					xc = randn(1,numP);
                    wave(m,:) = helmSolve(K,T,xc);
                end
            case 'trapezium'
                wave = zeros(numChn/2, numP+1);
                for m =1:numChn/2
                    % xc = pi*(2*rand(1,numP+1)-1);
					
					xc = randn(1,numP+1);
                    wave(m,:) = helmSolve(K,T,xc);
                end
        end

        if nargout >= 2
            wave_xy = zeros(numChn, size(wave,2));
            for m =1:numChn/2
                wave_xy(2*m-1:2*m,:) = [cos(wave(m,:)); sin(wave(m,:))];
            end
			if nargout == 3
				wave_xy_pl = zeros(numChn, size(wave,2));
				for m =1:numChn/2
					wave_xy_pl(2*m-1:2*m,:) =  control.pwr_levels{m} * [cos(wave(m,:)); sin(wave(m,:))];
				end
			end
        end

    case 'xy'
        switch control.integrator
            case 'rectangle'
                wave = zeros(numChn, numP);
                % for m =1:numChn
				
					% % maybe the error is senstive to fast-jumping shapes generetaed by sigmoid
                    % xc = 2*rand(1,numP)-1;
                    % xs = helmSolve(K,T,xc);
                    % wave(m,:) = sigmoid(xs);
					
					% % better error performance with linear and quadratic elements
					% % xc = randn(1,numP);
					% % wave(m,:) = helmSolve(K,T,xc);
                % end
				
				for m =1:numChn/2
				
					% better error performance with linear and quadratic elements
					xc = pi*(randn(1,numP));
					xs = helmSolve(K,T,xc);
					wave(2*m-1:2*m,:) = [cos(xs);sin(xs)];
					
                end

            case 'trapezium'
                wave = zeros(numChn, numP+1);
                % for m =1:numChn
                    % % xc = 2*rand(1,numP+1)-1;
                    % % xs = helmSolve(K,T,xc);
                    % % wave(m,:) = sigmoid(xs);
					
					% xc = randn(1,numP+1);
					% wave(m,:) = helmSolve(K,T,xc);
                % end
				
				for m =1:numChn/2
				
					% better error performance with linear and quadratic elements
					xc = pi*(randn(1,numP+1));
					xs = helmSolve(K,T,xc);
					wave(2*m-1:2*m,:) = [cos(xs);sin(xs)];
					
                end
        end

        if nargout >= 2
            wave_xy = wave;
			if nargout == 3
				wave_xy_pl = zeros(numChn, size(wave,2));
				for m =1:numChn/2
					wave_xy_pl(2*m-1:2*m,:) =  control.pwr_levels{m} * wave(2*m-1:2*m,:);
				end
			end
        end
end

end