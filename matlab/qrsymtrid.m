function [ Q, R ] = qrsymtrid( A )
% QR decomposition for symmetric tridiagonal matrices
% QRSYMTRID returns Q and R such that A=QR, Q is orthogonal (and upper
% Hessenberg), and R is upper triangular (of bandwidth 2, ie diagonal and
% two superdiagonals). R is returned as M x 3 matrix, with the bottom right
% corner empty (zero).

% reference: Golub, Van Loan; 3rd ed; 8.3.3
% $Id$

N = size( A, 1 );
if N ~= size( A, 2 ) || any(any(A ~= A'))
    error( 'input A must be symmetric' );
end;

R   = A;
Q   = eye( N );
CS  = zeros( N-1, 2 );

% remove 0 on subdiagonal by premultiplication with orthogonal Givens
% matrices Q1', Q2', ... -> obtain upper triangular R = QN-1'...Q2' Q1' T
% so T = QR with Q = Q1 Q2 Q3 ... QN-1
for k=1:(N-1)
    % now want to find Givens rotation in k, k+1 such that
    % R( k+1, k ) = 0. This needs to be applied to row k, k+1
    % GIVENS computes c and s such that [ c s; -s c ]' [ a b ]' = [ r 0 ]'
    [ c, s ]        = givens( R( k,k ), R( k+1, k ) );  % flop count: 6
    G               = [ c s; -s c ];
    K               = min( N, k+2 );                    % flop count: 1, or so...
    R( k:k+1, k:K ) = G' * R( k:k+1, k:K );             % flop count: 3*6
    R( k+1, k )     = 0;    % should be close to zero anyway...
    % accumulate Q
    Q(1:k+1, k:k+1) = Q(1:k+1, k:k+1) * G;         % flop count: (k+1)*6
end;

% flop count: without Q accumulation, around 25 N
% with Q accumulation: sum k=1:N of 30 + 6 k
% total flop count 3 N^2 + 33 N