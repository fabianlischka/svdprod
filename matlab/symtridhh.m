function A = symtridhh( A )
% SYMTRIDHH brings a symmetric real matrix A in similar tridiagonal form
% using similar householder transformations. 

if any(any(A ~= A'))
    error('Input must be symmetric!')
end
n = size( A, 1 );
for k = 1:(n-2)
    [ v, beta ] = house( A( (k+1):n, k ) );
    p = beta*A( (k+1):n, (k+1):n )*v;
    w = p - (beta*p'*v/2)*v;
    A(k+1,k) = norm( A( (k+1):n,k ) );  A( (k+2):n, k ) = 0;
    A(k,k+1) = A(k+1,k)              ;  A( k, (k+2):n ) = 0;
    % FIX: take advantage of symmetry here
    A( (k+1):n, (k+1):n ) = A( (k+1):n, (k+1):n ) - v*w' - w*v';
end;