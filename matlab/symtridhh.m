function [ A, Q ] = symtridhh( A )
% SYMTRIDHH brings a symmetric real matrix A in similar tridiagonal form
% using similar householder transformations, T = Q'AQ

% reference: Golub, Van Loan; 3rd ed; ch. 8.3.1 and 5.1.6.
% $Id$

if any(any(A ~= A'))
    error('Input must be symmetric!')
end
N = size( A, 1 );

for k = 1:(N-2)
    [ v, beta, mu ] = house( A( (k+1):N, k ) );                 % flops: 3(N-k)
    p = A( (k+1):N, (k+1):N )*(beta*v);                         % flops: 2(N-k)^2
    w = p - (beta/2*p'*v)*v;                                    % flops: 4(N-k)
    A( k, (k+2):N ) = 0              ;  A(k,k+1) = mu;          % mem:   N-k
    if nargout == 2
        % collect the v
        A( k+1:N, k )   = v          ;  A( k+1, k ) = beta; 
    else
        % else leave A symmetric, tridiagonal
        A( (k+2):N, k ) =            ;  A(k+1,k) = mu;          % mem:   N-k
    end;
    % FIX: take advantage of symmetry here
    A( (k+1):N, (k+1):N ) = A( (k+1):N, (k+1):N ) - v*w' - w*v';
    % flops: without taking advantage of sym:                   % flops: 3(N-k)^2
end;

if nargout == 2
    % accumulate Q
    Q = eye( N );
    for k = (N-2):-1:1
        v           = A( k+1:N, k );
        beta        = v( 1 );
        v( 1 )      = 1;
        Q(k+1:N, k+1:N) = Q( k+1:N, k+1:N ) - beta * v * v' * Q( k+1:N, k+1:N );
        % and fix the A
        A(k+1:N, k) = A(k, k+1:N )';
    end;
end;

% total flops: sum k=1:N of 5(N-k)^2 = sum k=1:N of 5N^2
% total flops: = 5/3 N^3 + lower order terms