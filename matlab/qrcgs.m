function [Q,R] = qrcgs( A )
% QR-decomposition, Classical Gram Schmidt
% (numerically unstable, Q might not be very orthogonal)

% reference: Golub, Van Loan; 3rd ed; 5.2.7
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
    % orthogonalisieren: compact version
    R(1:(k-1),k) = Q(:,1:(k-1))' * A(:,k);      % matrix product
    Q(:,k)       = A(:,k) - Q(:,1:(k-1)) * R(1:(k-1),k);
    % ie 4M*(k-1) flops
    
    % normalisieren
    R(k,k) = norm(Q(:,k),2);                    % flops: M mult, M-1 add, 1 sqrt
    % if vector empty, we need to create an orthogonal one
    if R(k,k) == 0
        % Q(k,k) = 1;
        % FIX FIX FIX
    else
        Q(:,k) = Q(:,k) / R(k,k);                   % flops: M mult
    end
    % another 3M flops
end

% flop count: sum k=1:N ( 3M + 4M*(k-1) ) = 3MN + 2MN^2
% for M>=N, the order is 2MN^2
% for M=N, 2 N^3