function [Q,R] = qrmgs( A )
% QR-decomposition, Modified Gram Schmidt (numerically stable, but Q might not be very orthogonal)
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
    % orthogonalisieren
    for n = 1:(k-1)
        R(n,k) = Q(:,n)' * Q(:,k);              % scalar product: M mult, M-1 add
        Q(:,k) = Q(:,k) - Q(:,n) * R(n,k);      % M mult, M add
    end
    % ie 4M*(k-1) flops
    
    % normalisieren
    R(k,k) = norm(Q(:,k),2);                    % M mult, M-1 add, 1 sqrt
    Q(:,k) = Q(:,k) / R(k,k);                   % M mult
    % another 3M flops
end

% flop count: sum k=1:N ( 3M + 4M*(k-1) ) = 3MN + 2MN^2
% for M>=N, the order is 2MN^2