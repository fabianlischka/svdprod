% QR-decomposition, Classical Gram Schmidt (numerically unstable)
% $Id$

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
K = min( M, N );    % min = steps

% important: we initialize Q with A
Q = A(:,1:K);       % Q will be M x K
R = zeros( K, N );  % R will be K x N

for k = 1:K
    % orthogonalisieren
    % iterative version
    if 0
    for n = 1:(k-1)
        R(n,k) = Q(:,n)' * A(:,k);              % scalar product: M mult, M-1 add
        Q(:,k) = Q(:,k) - Q(:,n) * R(n,k);      % M mult, M add
    end
    end
    
    % compact version
    R(1:(k-1),k) = Q(:,1:(k-1))' * A(:,k);      % matrix product
    Q(:,k)       = A(:,k) - Q(:,1:(k-1)) * R(1:(k-1),k);
    % ie 4M*(k-1) flops
    
    % normalisieren
    R(k,k) = norm(Q(:,k),2);                    % M mult, M-1 add, 1 sqrt
    Q(:,k) = Q(:,k) / R(k,k);                   % M mult
    % another 3M flops
end

% sum k=1:K ( 3M + 4M*(k-1) ) = 3MK + 2MK^2
% for M>=N, the order is 2MN^2