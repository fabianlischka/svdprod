function T = qrexstep( T );
% explicit QR iteration for sym matrices with Wilkinson shift
% note: we assume T to be symmetric, no check here

% reference: Golub, Van Loan; 3rd ed; 8.3.4
% $Id$

ident   = eye( size( T ) );                                                     
mu      = wilkinsonshift( T( end-1, end-1 ), T( end, end-1 ), T( end, end ) );                                                        
[Q,R]   = qrsymtrid( T - mu * ident );                  % flops: 3 N^2 + 33 N
% note: this can be optimized to take advantage of bandedness
T       = R*Q + mu * ident                              % flops: 
