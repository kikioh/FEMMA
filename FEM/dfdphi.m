% Generate an (N+1) x N sparse matrix A.
% where all non-zero elements are 1 and follow the given pattern,
% i.e., A(j,j) = 1, A(j+1,j) = 1 and the other elements in j-th column are zero
%
% mengjia.he@kit.edu, 2024.12.02


function A = dfdphi(N)

% Number of rows and columns
num_rows = N + 1;
num_cols = N;

% Row indices for non-zero elements (2 rows per column)
row_indices = [1:N; 2:N+1];

% Repeat the column indices for each block
col_indices = repmat(1:N, 2, 1);

% Convert row and column indices into vectors
row_indices = row_indices(:);
col_indices = col_indices(:);

% Create the sparse matrix
A = sparse(row_indices, col_indices, 1, num_rows, num_cols);

end
