% Generates operator for n-spin case. Syntax:
%
%   A = operatorCreate(operators,coherent,formalism)
%
% Parameters:
%
%   operators    - cell array, each element is a string, which gives the
%                state of invidival spin, can be 'Lx', 'Ly', 'Lz', 'L+', 'L-' or 'E'.
%                All the concerned spins should be specified, thosed
%                escaped spins should be specified with 'E'
%
%   coherent    - 'single', returns the sum of all single-spin operators,
%               or 'coherent', returns a product operator. The default choice 
%               is single.            
%
%   formalism   - 'hilb' for Hilbert space or 'liouv' for Liouville space.
%               The default choice is 'liouv'
%
% Outputs:
%
%   A           - sparse representation of spin operator (Hilbert space) or its
%               commutation superoperator (Liouville space).
%
% mengjia.he@kit.edu, 2025.02.28

function A = operatorCreate(operators,coherent,formalism)

% Default linear space
if ~exist('coherent','var'), coherent = 'single'; end

% Default linear space
if ~exist('formalism','var'), formalism = 'liouv'; end

% number of spins
numSpin = numel(operators);

% basis set
unit = sparse([1 0; 0 1]);
sigma_x = sparse([0 1; 1 0]/2);
sigma_y = sparse([0 -1i; 1i 0]/2);
sigma_z = sparse([1 0; 0 -1]/2);

% generate all single-spin operators
opt = cell(1,3*numSpin);
for n = 1 : numSpin
    Lx = 1; Ly = 1; Lz = 1;
    for k = 1 : numSpin
        if k == n
            Lx = kron(Lx, sigma_x); Ly = kron(Ly, sigma_y); Lz = kron(Lz, sigma_z);
        else
            Lx = kron(Lx, unit); Ly = kron(Ly, unit); Lz = kron(Lz, unit);
        end
    end
    opt{3*n-2} = Lx;    opt{3*n-1} = Ly;    opt{3*n} = Lz;  

end

% generate specified spin operator in Hilbert space
switch coherent

    case 'single'
        A = sparse(2^numSpin,2^numSpin);

        for n = 1:numSpin
            switch operators{n}

                case 'Lx', A = A + opt{3*n-2};
                case 'Ly', A = A + opt{3*n-1};
                case 'Lz', A = A + opt{3*n};    
                case 'L+', A = A + (opt{3*n-2} + 1i*opt{3*n-1});   
                case 'L-', A = A + (opt{3*n-2} - 1i*opt{3*n-1}); 
            end
        end

    case 'coherent'
        A = 1;
        for n = 1:numSpin
            switch operators{n}

                case 'Lx', A = A * opt{3*n-2};
                case 'Ly', A = A * opt{3*n-1};
                case 'Lz', A = A * opt{3*n};    
                case 'L+', A = A * (opt{3*n-2} + 1i*opt{3*n-1});   
                case 'L-', A = A * (opt{3*n-2} - 1i*opt{3*n-1}); 
            end
        end
end

% transform to Liouville space if required 
if strcmp(formalism,'liouv')

    % unit matrix
    unit_liouv = speye(2^(numSpin));
    A = kron(unit_liouv, A) - kron(transpose(A), unit_liouv);
end

A = full(A);
end