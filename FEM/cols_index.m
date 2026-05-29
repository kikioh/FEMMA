function cols = cols_index(N,nc,eleDof)


% Starting value of each block: [1, N+2, 2N+3, ...]
startVals = (0:nc-1)*N + (1:nc);

% Offsets inside each block: [0, 1, 2, ..., N-1]
offsets = 0:N-1;

% Use implicit expansion to generate all values
seq = startVals.' + offsets;   % nc×N matrix

% Convert to a single row vector
seq = seq.';                 % transpose to arrange by block

% element repeat to generate column index
cols = repelem(seq(:),eleDof);

end