function [ B, U, V ] = gksvdstep( B, U, V )
% Given an (upper) bidiagonal real B, NxN, with no zeros on the diagonal or
% superdiagonal, GKSVDSTEP (Golub-Kahan Singular Value Decomposition Step)
% returns B+  = U' B V, where U and V NxN are orthogonal matrices,
% and VV is the orthogonal matrix that would have been obtained by applying
% an implicit QR step with Wilkinson shift (QRIMSTEP) to T = B'B, 
% ie T+ = B+'B+ = V'B'U U'BV = V' B'B V = V'TV
% and we want T+ triangular, and V*e1 proportional to T-mu*I
%
% B is given in short format, ie. the diagonal is in the first col B(1:N,1), 
% the superdiagonal in the second col B(1:N-1, 2).

% reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.1
% $Id$

N = size( B, 1 );

% determine the shift mu.
% mu should be the eigenvalue (closer to the lower right) of
% ( B(N-1,1)^2 + B(N-2,2)^2    ,  B(N-1,1) * B(N-1,2)    )
% ( B(N-1,1) * B(N-1,2)        ,  B(N,1)^2 + B(N-1,2)^2  )
mu = wilkinsonshift( B(N-1,1)^2+B(N-2,2)^2,  B(N-1,1)*B(N-1,2),  B(N,1)^2+B(N-1,2)^2 );

x  = B( 1, 1 ) - mu;
y  = B( 1, 2 );

if nargin < 3
    V  = eye( N );
    if nargin < 2
        U  = V;
    end;
end;
    
% B+ = U' B V =   U(N-1)' ... U2' U1' B G1 V2 V3 ... V(N-1), with G1*e1 prop T-mu*I 
for k = 1:N-1
    % compute Givens rotation, around k,k+1 FROM RIGHT, G = [c s; -s c]
    [ c s ]     = givens( x, y );                       % flops: 5 + 1 sqrt    
    % rotate B: want B_(1/2) = B G, V+ = V G
    G           = [ c s; -s c ];
    K           = min( N, k+2 );                        % flops: 1, or so...
    % update V
    V( 1:K, k:k+1 ) = V( 1:K, k:k+1 ) * G;              % flops: 6 * K
    % for B in ordinary format, we want:
    % B( k-1:K, k:k+1 ) = B( k-1:K, k:k+1 ) * G;

    if k > 1
        % for debugging we can compute the zero'ed out element:
        zerod   = s * B(k-1,2) + c * bulge
        B( k-1, 2 ) =  c * B( k-1, 2 ) - s * bulge;     % flops: 3
    end;
    Bk          = B( k, 1 );
    bulge       = -s * B( k+1, 1 );                     % flops: 1
    B( k, 1 )   =  c * Bk - s * B( k, 2 );              % flops: 3
    B( k, 2 )   =  s * Bk + c * B( k, 2 );              % flops: 3
    B( k+1, 1 ) =  c * B( k+1, 1 );                     % flops: 1

    % now, bulge is at k+1, k
    x           = B( k, 1 );
    y           = bulge;       % which we want to eliminate...
    
    % compute Givens rotation, around k,k+1 FROM LEFT, G = [c s; -s c]
    [ c s ] = givens( x, y );                           % flops: 5 + 1 sqrt    
    % rotate B: want B_+ = G' B_(1/2), U+ = U G
    G               = [ c s; -s c ];
    % update U
    U( 1:K, k:k+1 ) = U( 1:K, k:k+1 ) * G;              % flops: 6 * K
    % for B in ordinary format, we want:
    % B = G' * B;

    zerod       =  s * B( k, 1 ) + c * bulge
    B( k, 1 )   =  c * B( k, 1 ) - s * bulge;
    Bk2         = B( k, 2 );
    B( k, 2 )   =  c * Bk2       - s * B( k+1, 1 );
    bulge       = -s * B( k+1, 2 );
    B( k+1, 1 ) =  s * Bk2       + c * B( k+1, 1 );
    B( k+1, 2 ) =  c * B( k+1, 2 );
    % now, bulge is at k, k+2

    x           = B( k, 2 );
    y           = bulge;       % which we want to eliminate...
end;
B( N, 2 ) = 0;