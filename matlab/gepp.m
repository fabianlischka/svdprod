% Gaussian Elimination, partial pivoting

M = size( A, 1 )    % rows
N = size( A, 2 )    % cols
K = min( M, N )     % min = steps

U = A;              % U will be K x N
L = eye( M, K );    % L will be M x K

% iterative version: k=steps m=rows n=cols
if 0
for k = 1:(K-1)
    % choose pivot
    [ Y, Index ] = max(abs(U(k,:)));
    % pivot is U(k,Index)
    % iterate through rows below
    for m = (k+1):M
        mult = U(m, Index) / U(k, Index);
        L(m, k) = mult;
        for n = 1:N
            U(m, n) = U(m, n) - mult * U(k, n);
        end
    end
end
end


% vectorized version: k m n
if 1
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