% TESTQR runs a few QR decompositions for a test matrix, 
% and reports the orthogonality (norm of Q'Q-I), and norm of A-QR

% $Id$

qralgos = { 'qr', 'qrhh', 'qrcgs', 'qrmgs' };

N = 20;
disp( sprintf( '\nTesting various QR algorithms, reporting orthogonality (ie norm of Q''Q-I)\n and norm of A-QR (smaller=better)')); 
for Typ = 1:6
	TestMat = gentestmat( N, Typ );
	
	for k = 1:size(qralgos,2)
        [Q,R]   = feval( qralgos{ k }, TestMat );
        Ortho   = norm( Q'*Q - eye( N ) );
        ErrNorm = norm( TestMat-Q*R );
        disp( sprintf( 'Algo %6s: Ortho %12g, ErrNorm %12g', qralgos{k}, Ortho, ErrNorm ) );
	end;
end;