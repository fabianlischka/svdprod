function [ B, U, V ] = bidighh( A );
% BIDIGHH computes bidiagonal B=U'AV using householder reflections
% suppose A MxN, M >= N, A real. This routine finds orthogonal U MxM, V NxN 
% such that B = U'AV is bidiagonal (the diagonal and one superdiagonal).
% The diagonal is stored in the first col of B, the superdiagonal in the
% second col of B, and B(end,2) == 0

% reference: Golub, Van Loan; 3rd ed; 5.4.3
% $Id$

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
if M < N
    error( 'A must have M >= N' );
end;

B       = A;              % B will be M x N
betas   = zeros( N, 1 );
gammas  = zeros( N-2, 1 );

for k = 1:N
    % determine HH reflector for column B(k:M,k)
    [ v, beta, mu ] = house( B(k:M,k) );            % flops: 3(M-k) 
    betas( k )      = beta;
    % B(k,k)          = mu;                           % = norm(x)
    B(k:M,k:N)      = B(k:M,k:N) - beta * v * v' * B(k:M,k:N);
                                                    % flops: 4*(N-k)*(M-k)
    B((k+1):M,k)    = v(2:end);
    % or, if you don't want to collect it, B((k+1):M,k) = 0;
    % B((k+1):M,k)    = 0;
    if k < N - 1
		[ v, beta, mu ] = house( B(k,(k+1):N)' );   % flops: 3(N-k)
		gammas( k )     = beta;
		% B(k,k+1)        = mu;                        % = norm(x)
        % note that we need to multiply from the right here:
		B(k:M,(k+1):N)  = B(k:M,(k+1):N) - beta * B(k:M,(k+1):N) * v * v';
                                                    % flops: 4*(N-k)*(M-k)
		B(k,(k+2):N)    = v(2:end)';
        % or, if you don't want to collect it, B(k,(k+2):N) = 0;
        % B(k,(k+2):N)    = 0;
    end
end;

% total flops (without computing U, V)
% sum k = 1:N of 8(N-k)(M-K) + 6(M-k)
% = 4 MN^2 - 4/3 N^3 = 4 N^2 (M-N/3), or, 
% for M=N, 8/3 N^3

if nargout > 1
    % compute U and V
    v = zeros( M, 1 );
    U = eye( M );
    for k = N:-1:1
        % HH reflector now: I - beta * v * v'. Note: v(1) = 1
        v(k)        = 1;
        v(k+1:M)    = B((k+1):M,k);
        U(k:M,k:M)  = U(k:M,k:M) - betas( k ) * v(k:M) * v(k:M)' * U(k:M,k:M);
        % flop count: 4*(M-k)^2
        B(k+1:M,k)  = 0;
    end;
    V = eye( N );
    for k = N-1:-1:2
        v(k)        = 1;
        v(k+1:N)    = B(k-1,(k+1):N)';
        V(k:N,k:N)  = V(k:N,k:N) - gammas( k-1 ) * v(k:N) * v(k:N)' * V(k:N,k:N);
        % flop count: 4*(N-k)^2
        B(k-1,k+1:N)= 0;
    end;
end;
B = [ diag(B) [ diag(B,1); 0 ] ];

% flop count for V: 4/3 N^3
% flop count for U: 4( NM^2 - MN^2 + N^3/3 )
% total flops (with computing U,V)
% 4 NM^2 + 4/3 N^3, or for M==N, 16/3 N^3 = 5.3 N^3