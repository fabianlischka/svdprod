function A = symtridhh( A )
% SYMTRIDHH brings a symmetric real matrix A in similar tridiagonal form
% using similar householder transformations. 

% $Header$

if any(any(A ~= A'))
    error('Input must be symmetric!')
end
N = size( A, 1 );
for k = 1:(N-2)
    [ v, beta, mu ] = house( A( (k+1):N, k ) );                 % flops: 3(N-k)
    p = A( (k+1):N, (k+1):N )*(beta*v);                         % flops: 2(N-k)^2
    w = p - (beta/2*p'*v)*v;                                    % flops: 4(N-k)
    A(k+1,k) = mu                    ;  A( (k+2):N, k ) = 0;    % mem:   N-k
    A(k,k+1) = A(k+1,k)              ;  A( k, (k+2):N ) = 0;    % mem:   N-k
    % FIX: take advantage of symmetry here
    A( (k+1):N, (k+1):N ) = A( (k+1):N, (k+1):N ) - v*w' - w*v';
    % flops: without taking advantage of sym:                   % flops: 3(N-k)^2
end;

% total flops: sum k=1:N of 5(N-k)^2 = sum k=1:N of 5N^2
% total flops: = 5/3 N^3 + lower order terms