% Gaussian Elimination, partial pivoting
% $Id$

M = size( A, 1 )    % rows
N = size( A, 2 )    % cols
K = min( M, N )     % min = steps

U = A;              % U will be K x N
L = zeros( M, K );  % L will be M x K
P = 1:K;            % this contains the permutation: row k has been permuted with P(k)

% iterative version: k=steps m=rows n=cols
if 1
for k = 1:(K-1)
    P
    L
    U
    % choose pivot
    [ Y, IndirectIndex ] = max(abs(U(P(k:K),k)));
    Index = P( k-1+IndirectIndex )
    % pivot is U(Index,k)
    temp        = P(Index);
    P(Index)    = P(k);
    P(k)        = temp;
    % note: now P(1:(k-1)) contains rows that have been done, P(k) == Index
    % contains the row we do now, P((k+1):K) contains untouched rows (which
    % need to be touched)
    L( Index, k ) = 1;
    % iterate through other rows: these are in P((k+1):K)
    for m = P((k+1):K)
        mult = U(m, k) / U( Index, k )   % note: Y == abs( U( Index, k ) ), ie Y is abs(pivot), not the pivot value itself...
        L(m, k) = mult;
        for n = 1:N     % need not be 1:N - need to think
            U(m, n) = U(m, n) - mult * U(Index, n);
        end
    end
    P
    L
    U
end
end


% vectorized version: k m n
if 0
for k = 1:(K-1)
    % choose pivot
    [ Y, Index ] = max(abs(U(k,:)));
    % pivot is U(k,Index)
    mult = U((k+1):M, Index) / U(k, Index);
    L((k+1):M, k ) = mult;
    for m = (k+1):M
        U(m, 1:N) = U(m, 1:N) - mult(m-k) * U(k, 1:N);
    end
end
end