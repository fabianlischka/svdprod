function [ B, U, V ] = bidigprod( A )
% BIDIGHH computes bidiagonal B=U' AAA V using householder reflections.
% suppose A NxNxK, real. Without ever computing the product 
% AAA = A(:,:,K)*...*A(:,:,2)*A(:,:,1), this routine finds 
% orthogonal U NxN, V NxN such that B = U' AAA V is bidiagonal 
% (the diagonal and one superdiagonal).
% B is returned in short format, ie the diagonal
% is stored in the first col of B, the superdiagonal in the
% second col of B, and B(end,2) == 0

% reference: Golub, Solna, van Dooren: Computing the SVD
% of a General Matrix Product/Quotient, 
% SIAM J. MATRIX ANAL. APPL., Vol. 22, No. 1, pp 1-19
% $Id$

N   = size( A, 1 );
if size( A, 2 ) ~= N
    error( 'matrices in A must be square (ie size(A,1)=size(A,2))' )
end;
    
K   = size( A, 3 );
V   = eye( N );
A(:,:,K+1)  = V;    % this is U

for t=1:N-1
    for k=1:K
        [ v, beta, mu ] = house( A( t:N, t, k ) );
        % I-beta vv' is Q^t_k, applied to t:N
        % multiply Ak from left...
        A( t:N, t:N, k )    = A( t:N, t:N, k   ) - beta * v * v' * A( t:N, t:N, k );
        % and Ak+1 from right - here we need to multiply all columns,
        % since we don't have the zeros needed in here
        A(  : , t:N, k+1 )  = A(  : , t:N, k+1 ) - A(  : , t:N, k+1 ) * v * v' * beta;
    end;

    if t < N-1
        % now, determine the t-th row of the product
        row = A( t, t:N, K );
        for k = K-1:-1:1
            row = row * A( t:N, t:N, k );
        end;
        % and determine householder to eliminate it!
        % (the diagonal element remains untouched (and is thrown away))
        [ v, beta, mu ] = house( row( 1, 2:end )' );    
        % multiply A( 1 ) and V from the right
        % need to multiply full columns (all rows), since might be filled
        A( :, t+1:N, 1 )    = A( :, t+1:N, 1 ) - A( :, t+1:N, 1 ) * v * v' * beta;
        V( :, t+1:N )       = V( :, t+1:N )    - V( :, t+1:N )    * v * v' * beta;
    end;
end;

% now B = A(:,:,K+1)'*AA*V = A(:,:,K)*...*A(:,:,2)*A(:,:,1)
q   = ones(N,1);
e   = zeros(N-1,1);

for k=1:K
    d   = diag( A(:,:,k) );
    e   = e .* d(1:end-1) + q(2:N) .* diag( A(:,:,k), 1 );
    q   = q .* d;
end;

U = A(:,:,K+1);
B = [ q [e; 0] ];
