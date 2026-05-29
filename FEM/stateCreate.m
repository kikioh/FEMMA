% generate the spin state under Zeeman basis. Syntax:
%
%           rho = stateCreate(states,coherent,formalism)
%
% Parameters:
%
%   states     - cell array, each element is a string, which gives the
%                state of invidival spin, can be 'Lz', 'Lx', 'Ly', 'L+', 'L-' or 'E'.
%                All the concerned spins should be specified.
%
%   coherent    - 'single', returns the sum of all single-spin matrices/vectors
%               - 'coherent', returns a product state density matrix/vector. 
%               The default choice is single.            
%
%   formalism   - 'hilb' for Hilbert space or 'liouv' for Liouville space.
%                  The default choice is 'liouv'
%
% Output
%
%   rho         - a Hilbert space density matrix or a Liouville space state vector
%
% mengjia.he@kit.edu, 2025.02.28

function rho = stateCreate(states,coherent,formalism)

% Default linear space
if ~exist('coherent','var'), coherent = 'single'; end

% Default linear space
if ~exist('formalism','var'), formalism = 'liouv'; end

% generate spin states in Hilbert space
rho = operatorCreate(states,coherent,'hilb');

% transform to Liouville space if required
if strcmp(formalism,'liouv'), rho = rho(:); end

% normalize state
rho = rho/norm(rho,'fro');

end