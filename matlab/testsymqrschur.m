% $Id$

N = 60;

disp( 'We compute eigenvalues of a test matrix using eig(), and then using the schur decomposition' );
disp( 'computed by symqrschur, and display the norm of the difference of the vector of eigenvalues' );
disp( 'relative to the norm of the vector of eigenvalues itself. Also the norm of A-QDQ'', and Q''Q-I.' );
disp( 'All errors should be around 1e-15.' );
disp( sprintf( '\nTesting SYMQRSCHUR with various (symmetrized) %g x %g matrices', N, N )); 
for Typ = 1:14
	TestMat = gentestmat( Typ, N, N, 0 );
    TestMat = TestMat + TestMat';
    Ev      = eig( TestMat );
    [ D,Q ] = symqrschur( TestMat );
    ErrEv   = sort( D ) - sort( Ev );
    Res     = norm( TestMat - Q*diag(D)*Q' );
    RelRes  = Res / norm( TestMat );
    Orth    = norm( Q'*Q - eye( N ) );
    % disp( sprintf( 'eval 1norm: %12g, infnorm: %12g; Res: %12g', norm( ErrEv, 1 ), norm( ErrEv, 'inf' ), Res ) );
    disp( sprintf( 'Type %2g: eval 1norm: %12g, infnorm: %12g; Res: %12g, Orth: %12g', Typ, norm(ErrEv,1)/norm(Ev,1), norm(ErrEv,'inf')/norm(Ev,'inf'), RelRes, Orth ) );
end;