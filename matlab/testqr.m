% TESTQR runs a few QR decompositions for a test matrix, 
% and reports the orthogonality (norm of Q'Q-I), and norm of A-QR
% $Id$

qralgos = { 'qr', 'qrhh', 'qrcgs', 'qrmgs' };

N=60;

% matrix with very different singular values: 1/2 .. 1/(2^N)
[U,X]=qr(rand(N));
[V,X]=qr(rand(N));
%TestMat = U* diag(2.^(-1:-1:-N)) *V;

% matrix with numbers of very different size
%TestMat = exp(16*rand(N)-20).*(rand(N)-.5);

% rank 1 matrix with slight disturbance
%TestMat = 5 * repmat( randn(N,1), 1, N ) + 1e-15*rand(N);   

TestMat = U* ( diag(1e-10*ones(N,1)) + diag(ones(N-1,1),1) ) * U';

for k = 1:size(qralgos,2)
    [Q,R] = feval( qralgos{k}, TestMat );
    Ortho   = norm( Q'*Q - eye( N ) );
    ErrNorm = norm( TestMat-Q*R );
    disp( sprintf( 'Algo %6s: Ortho %12g, Norm %12g', qralgos{k}, Ortho, ErrNorm ) );
end;