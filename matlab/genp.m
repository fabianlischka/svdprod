% Gaussian Elimination, no pivoting
% $Id$

M = size( A, 1 )    % rows
N = size( A, 2 )    % cols
K = min( M, N )     % min = steps

U = A;              % U will be K x N
L = eye( M, K );    % L will be M x K

% iterative version: k m n
if 1
for k = 1:K         % steps
    % pivot is A(k,k)
    for m = (k+1):M     % iterate through rows below, determine multiplier, apply to row
        mult   = U(m, k) / U(k, k);     % 1 flop, M-k-1 times
        L(m,k) = mult;
        for n = k:N     % iterate through cols to right, could also be 1:N, same result, slower
            U(m, n) = U(m, n) - mult * U(k, n); % 2 flops, N-k times
        end
    end
end
end
% Complexity: N^2 * ( M - N + 2/3 N) = N^2 * (M - 1/3 N) [ = 2/3 N^3 for M=N ] 
% 2/3 N^3
U=U(1:K,:);

% vectorized version: k m n
if 0
for k = 1:(K-1)
    % pivot is A(k,k)
    mult = U((k+1):M, k) / U(k, k);
    L((k+1):M, k ) = mult;
    for m = (k+1):M
        U(m, k:N) = U(m, k:N) - mult(m-k) * U(k, k:N);
    end
end
end
