function [Q,R] = qrcgs( A )
% QR-decomposition, Classical Gram Schmidt (numerically unstable)
% $Id$

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
if M < N
    err( 'A must have number rows >= number columns' );
end;

% important: we initialize Q with A
Q = A;              % Q will be M x N
R = zeros( N, N );  % R will be N x N

for k = 1:N
    % orthogonalisieren: iterative version, old
    % for n = 1:(k-1)
    %     R(n,k) = Q(:,n)' * A(:,k);              % scalar product: M mult, M-1 add
    %     Q(:,k) = Q(:,k) - Q(:,n) * R(n,k);      % M mult, M add
    % end
    
    % orthogonalisieren: compact version
    R(1:(k-1),k) = Q(:,1:(k-1))' * A(:,k);      % matrix product
    Q(:,k)       = A(:,k) - Q(:,1:(k-1)) * R(1:(k-1),k);
    % ie 4M*(k-1) flops
    
    % normalisieren
    R(k,k) = norm(Q(:,k),2);                    % M mult, M-1 add, 1 sqrt
    Q(:,k) = Q(:,k) / R(k,k);                   % M mult
    % another 3M flops
end

% flop count: sum k=1:N ( 3M + 4M*(k-1) ) = 3MN + 2MN^2
% for M>=N, the order is 2MN^2