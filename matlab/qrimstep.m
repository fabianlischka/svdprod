% symmetric implicit QR step with Wilkinson shift
% given a symmetric tridiagonal T, routine computes T+ = Q'TQ, 
% where QR = T - mu I  is a QR decomposition of T shifted by mu

% Q = G1 G2 .. G(N-1), where we determine Givens rotations ..
% .. G1 such that e1*G1 parallel to first col of T-mu I
% .. G2, ..., G(N-1) such that bulge above is chased down and out

% reference: Golub, Van Loan; 3rd ed; ch. 8.3.5
% $Id$

    	
    TS = stridl2s( TL );
    N = size( TL, 1 );
    % note: T in short format
    
    d       = ( TS( N-1, 1 ) - TS( N, 1 ) ) / 2;    % flops: 2
    % mu      = TS( N, 1 ) + d - sign( d ) * sqrt( d^2 + TS(N-1, 2)^2 );
    % note that ( d - sign(d)*sqrt(d^2+b^2) ) * ( d + sign(d)*sqrt(d^2+b^2) ) =
    %           = d^2 - (d^2 + b^2) = - b^2, so we can avoid cancelation
    %           between d and the square root, particularly when b^2 (which
    %           is the last superdiagonal element) is small - which we hope!
    mu      = TS( N, 1 ) - ( TS(N-1, 2)^2 / ( d + sign(d)*sqrt(d^2 + TS(N-1, 2)^2) ) );

    %   note: Givens(x) gives us G with G'x=scalar * e1,
	%   so GG'x = x = scalar * G * e1, so we can chose x=(T-mu I)e1, 
	%   and use the normal Givens routine!
    x       = TS(1, 1) - mu;
    y       = TS(1, 2);
    
    % note: later, we would want to pass a matrix in to be updated
    % (instead of starting from the identity, and then multiply...)
    Q       = eye( N );
    
    for k = 1:(N-1)
        % compute Givens rotation, around k,k+1, G = [c s; -s c]
        [ c s ] = givens( x, y );                           % flops: 5 + 1 sqrt
        % rotate T: want T+ = G' T G, and Q+ = Q G
        G               = [ c s; -s c ];
        K               = min( N, k+2 );                    % flops: 1, or so...
        Q( 1:K, k:k+1 ) = Q( 1:K, k:k+1 ) * G;              % flops: 6 * K
        % for long T, we want:
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
        TS( k, 1 )   = c^2 * Tk1 + s^2 * TS( k+1, 1 ) - 2*c*s*Tk2;
        TS( k, 2 )   = c*s*( Tk1 - TS( k+1, 1 ) ) + ( c^2 - s^2 ) * Tk2;
        TS( k+1, 1 ) = s^2 * Tk1 + c^2 * TS( k+1, 1 ) + 2*c*s*Tk2;
        TS( k+1, 2 ) = c * TS( k+1, 2 );
                                            % last 4 lines: flops: +- 27        
        x            = TS( k, 2 );
        y            = bulge;       % which we want to eliminate...
    end;
    TS( N, 2 ) = 0;

    % flop count: around 40 N for new TS, around 3 N^2 for accumulating Q