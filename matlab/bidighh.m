% bidiagonalization using householder
% suppose A MxN, M >= N, A real. We will find orthogonal U MxM, V NxN 
% such that R = U'AV is bidiagonal (the diagonal and one superdiagonal)
% note: U, V not actually computed here

% reference: Golub, Van Loan; 3rd ed; 5.4.3
% $Id$

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
if M < N
    err( 'A must have M >= N' );
end;

R       = A;              % R will be M x N
betas   = zeros( N-1, 1 );
gammas  = zeros( N-2, 1 );

for k = 1:N
    % determine HH reflector for column A(k:N,k)
    [ v, beta, mu ] = house( R(k:M,k) );            % flops: 3(M-k) 
    betas( k )      = beta;
    R(k,k)          = mu;                           % = norm(x)
    R(k:M,(k+1):N)  = R(k:M,(k+1):N) - beta * v * v' * R(k:M,(k+1):N);
                                                    % flops: 4*(N-k)*(M-k)
    % R((k+1):M,k)    = v(2:end);
    % or, if you don't want to collect it, R((k+1):M,k) = 0;
    R((k+1):M,k)    = 0;
    if k < N - 1
		[ v, beta, mu ] = house( R(k,(k+1):N)' );   % flops: 3(N-k)
		gammas( k )     = beta;
		R(k,k+1)        = mu;                        % = norm(x)
        % note that we need to multiply from the right here:
		R((k+1):M,(k+1):N)  = R((k+1):M,(k+1):N) - beta * R((k+1):M,(k+1):N) * v * v';
                                                    % flops: 4*(N-k)*(M-k)
		% R(k,(k+2):N)    = v(2:end)';
        % or, if you don't want to collect it, R(k,(k+2):N) = 0;
        R(k,(k+2):N)    = 0;
    end
end;

% total flops (without computing U, V)
% sum k = 1:N of 8(N-k)(M-K) + 6(M-k)
% = 4 MN^2 - 4/3 N^3 = 4 N^2 (M-N/3), or, 
% for M=N, 8/3 N^3

% if desired, compute Q, and zero out R
%for k = (N-1):-1:1
    % HH reflector now: I - beta * v * v'. Note: v(1) = 1
%    v(k)=1;
%    v((k+1):N) = R((k+1):N,k);
%    Q(k:N,k:N) = Q(k:N,k:N) - betas( k ) * v(k:N) * transpose(v(k:N)) * Q(k:N,k:N);
%    R((k+1):N,k) = 0;
%end;