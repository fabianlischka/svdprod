function [ D, U, V ] = gksvd( A, tol )
% GKSVD computes the SVD of A
% Given A MxN, M>=N, and tol, this computes U MxM, V NxN orthogonal and DD diagonal 
% such that A=U(DD+E)V', where two norm E is around unit roundoff x two-norm A
% D is returned as a Nx1 column vector, DD would be = [ diag(D); zeros( M-N, N ) ];

% reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.2
% $Id$

if nargin == 1
    tol = 1e-14;
end;

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
if M < N
    error( 'A must have M >= N' );
end;

% BIDIGHH computes bidiagonal B=U'AV using householder reflections
[ B, U, V ] = bidighh( A );
% GKSVDSTEPS computes diagonal D and updates U, V such that D=U'AV
[ D, U, V ] = gksvdsteps( B, U, V, tol );
