function [ TS, Q ] = qrimstep( TS, Q ) 
% symmetric implicit QR step with Wilkinson shift
% given a symmetric tridiagonal T, routine computes T+ = Q'TQ, 
% where QR = T - mu I  is a QR decomposition of T shifted by mu

% Q = G1 G2 .. G(N-1), where we determine Givens rotations ..
% .. G1 such that e1*G1 parallel to first col of T-mu I
% .. G2, ..., G(N-1) such that bulge above is chased down and out

% reference: Golub, Van Loan; 3rd ed; ch. 8.3.5
% $Id$

N = size( TS, 1 );
% note: T in short format

% determine Wilkinson shift
mu = wilkinsonshift( TS( N-1, 1 ), TS( N-1, 2 ), TS( N, 1 ) );
    % flop count: 15

%   note: Givens(x) gives us G with G'x=scalar * e1,
%   so GG'x = x = scalar * G * e1, so we can chose x=(T-mu I)e1, 
%   and use the normal Givens routine!
x       = TS(1, 1) - mu;
y       = TS(1, 2);

% note: later, we would want to pass a matrix in to be updated
% (instead of starting from the identity, and then multiply...)
% note: in flop count, we assume Q has M rows
if nargin < 2
    Q   = eye( N );
end;

for k = 1:(N-1)
    % compute Givens rotation, around k,k+1, G = [c s; -s c]
    [ c s ]         = givens( x, y );                   % flops: 5 + 1 sqrt
    % rotate T: want T+ = G' T G, and Q+ = Q G
    G               = [ c s; -s c ];
    K               = min( N, k+2 );                    % flops: 1, or so...
    % update Q
    Q(:, k:k+1 )    = Q(:, k:k+1 ) * G;                 % flops: 6 * M
    % for T in ordinary format, we want:
    % T( k:k+1, k-1:K ) = G' * T( k:k+1, k-1:K );
    % T( k-1:K, k:k+1 ) = T( k-1:K, k:k+1 ) * G;
    
    if k > 1
        % for debugging we can compute the zero'ed out element:
        % zerod = s * TS(k-1,2) + c * bulge
        TS( k-1, 2 ) =  c * TS( k-1, 2 ) - s * bulge;   % flops: 3
    end;
    Tk1          = TS( k, 1 );
    Tk2          = TS( k, 2 );
    bulge        = -s * TS( k+1, 2 );                   % flops: 1
    c2           = c^2;  s2 = s^2; sc2 = 2*c*s;
    TS( k, 1 )   = c2  * Tk1 + s2 * TS( k+1, 1 ) - sc2*Tk2;
    TS( k, 2 )   = c*s*( Tk1 - TS( k+1, 1 ) ) + ( c2 - s2 ) * Tk2;
    TS( k+1, 1 ) = s2  * Tk1 + c2 * TS( k+1, 1 ) + sc2*Tk2;
    TS( k+1, 2 ) = c   * TS( k+1, 2 );
    % last 4 lines: flops: +- 27, could be optimized (by storing c^2 etc.)
    % note: should seek to avoid cancellation above
    x            = TS( k, 2 );
    y            = bulge;       % which we want to eliminate...
end;
TS( N, 2 ) = 0;

% flop count: around 32 N for new TS, around 6 MN for accumulating Q