% Builds auxiliary matrices for the calculation of the directional 
% derivatives of the trapezium product quadrature propagator:
%
%           expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt)
% 
% with respect to control coefficients in the evolution generators
% HL and HR on the left and the right edge of the interval. The de-
% rivatives are calculated using Eq 16 of Goodwin and Kuprov:
%
%                 https://doi.org/10.1063/1.4928978
%
% Syntax:
%
%        [G_fwd_l,G_bwd_l,G_fwd_r,G_bwd_r]=generator_trap(drifts,controls,dt,cL,cR,k,h)
%
% Parameters:
%
%     drifts - a cell array of two matrices containing drift 
%              generators at the left (first element) and the
%              right (second element) edge of the interval
%
%   controls - a cell array of control generators
%         
%         dt - interval duration, seconds
%
%         cL - control generator coefficients at the left 
%              edge of the interval
%
%         cR - control generator coefficients at the right
%              edge of the interval
%
%          k - the number of the generator inside controls 
%              array that the differentiation refers to
%
%		   h - Finite difference step
%
% Outputs:
%
%     G_fwd_l 	- left-hand forward generator 
%
%     G_bwd_l 	- left-hand backward generator 
%
%     G_fwd_r 	- right-hand forward generator 
%
%     G_bwd_r 	- right-hand backward generator 
%
%
% u.rasulov@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% This function is a modified version of the original from the Spinach 
% library (MIT License). Original source: 
% <https://spindynamics.org/wiki/index.php?title=aux_mat.m>
%
% Modified by mengjia.he@kit.edu, 2025.11.07

function [G_fwd_l,G_bwd_l,G_fwd_r,G_bwd_r]=generator_trap(drifts,controls,dt,cL,cR,k,h)

% Check consistency
grumble(drifts,controls,dt,cL,cR,k);

% Build left and right generators
G{1}=drifts{1}; G{2}=drifts{2};
for n=1:numel(controls)
    G{1}=G{1}+cL(n)*controls{n};
    G{2}=G{2}+cR(n)*controls{n};
end

% Build the overall generator
G_fwd_l = G{1} + h * controls{k};
G_fwd_l = (G_fwd_l+G{2})/2+1i*dt*(1/12)*(G_fwd_l*G{2}-G{2}*G_fwd_l);

G_bwd_l = G{1} - h * controls{k};
G_bwd_l = (G_bwd_l+G{2})/2+1i*dt*(1/12)*(G_bwd_l*G{2}-G{2}*G_bwd_l);

G_fwd_r = G{2} + h * controls{k};
G_fwd_r = (G{1}+G_fwd_r)/2+1i*dt*(1/12)*(G{1}*G_fwd_r-G_fwd_r*G{1});

G_bwd_r = G{2} - h * controls{k};
G_bwd_r = (G{1}+G_bwd_r)/2+1i*dt*(1/12)*(G{1}*G_bwd_r-G_bwd_r*G{1});


end

% Consistency enforcement
function grumble(drifts,controls,dt,cL,cR,k)
if ~iscell(drifts)
    error('drifts must be a cell array of matrices.');
end
for n=1:numel(drifts)     
    if (~isnumeric(drifts{n}))||(size(drifts{n},1)~=size(drifts{n},2))
        error('all elements of drifts cell array must be square matrices.');
    end
end
if ~iscell(controls)
    error('controls must be a cell array of square matrices.');
end
for n=1:numel(controls)
    if (~isnumeric(controls{n}))||...
       (size(controls{n},1)~=size(controls{n},2))||...
       (size(controls{n},1)~=size(drifts{1},1))
        error('control generators must have the same size as drift generators.');
    end
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(cL))||(~isreal(cL))||(~isnumeric(cR))||(~isreal(cR))
    error('cL and cR must be arrays of real numbers.');
end
if (~isnumeric(k))||(~isreal(k))||(~isscalar(k))||(k<=0)||(mod(k,1)~=0)
    error('k must be a positive real integer.');
end
end

% Rutherford to Eddington (who wondered at table whether we should ever
% come to know electrons as more than mental concepts): "Not exist? Not
% exist?! Why, I can see the little buggers as plain as I can see that 
% spoon in front of me!"

