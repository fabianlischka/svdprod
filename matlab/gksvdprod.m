function [ D, U, V ] = gksvdprod( A, tol )
% GKSVDPROD computes the SVD of a product of square matrices
% Given A NxNxK, for AAA := A(:,:,K)*...*A(:,:,1) this computes
% U, V NxN orthogonal, and D Nx1 such that AAA=U diag(D) V'
% Note: AAA NxN is never explicitly computed!

% reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.2
% reference: Golub, Solna, van Dooren: Computing the SVD
% of a General Matrix Product/Quotient, 
% SIAM J. MATRIX ANAL. APPL., Vol. 22, No. 1, pp 1-19
% $Id$

if nargin == 1
    tol = 1e-14;
end;

% BIDIGPROD computes bidiagonal B=U' AAA V using householder reflections
[ B, U, V ] = bidigprod( A );
% GKSVDSTEPS computes diagonal D and updates U, V such that D=U'AV
[ D, U, V ] = gksvdsteps( B, U, V, tol );
