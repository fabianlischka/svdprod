function T = qrexstep( T );
% explicit QR iteration for sym matrices with Wilkinson shift

% reference: Golub, Van Loan; 3rd ed; 8.3.4
% $Id$


ident   = eye( size( T ) );
% T       = symtridhh( A );
% mu  = A( end, end )
d       = ( T( end-1, end-1 ) - T( end, end ) ) / 2;    % flops: 2
mu      = T( end, end ) + d - sign( d ) * sqrt( d^2 + T(end, end-1)^2 )
                                                        % flops: ~8
[Q,R]   = qrsymtrid( T - mu * ident );                  % flops: 3 N^2 + 33 N
% note: this can be optimized to take advantage of bandedness
T       = R*Q + mu * ident                              % flops: 
