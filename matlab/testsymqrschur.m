% $Id$

N = 60;

disp( 'We compute eigenvalues of a test matrix using eig(), and then using the schur decomposition' );
disp( 'computed by symqrschur, and display the norm of the difference of the vector of eigenvalues' );
disp( 'and also relative to the norm of the vector of eigenvalues itself. The latter should always' );
disp( 'be around 1e-15.' );
disp( sprintf( '\nTesting SYMQRSCHUR with various (symmetrized) %g x %g matrices', N, N )); 
for Typ = 1:14
	TestMat = gentestmat( Typ, N, N, 0 );
    TestMat = TestMat + TestMat';
    Ev      = eig( TestMat );
    ErrEv   = sort(symqrschur( TestMat ) ) - sort( Ev );
    disp( sprintf( 'for type %g', Typ ) );
    disp( sprintf( 'err norms: 1norm: %12g,  2norm: %12g,  infnorm: %12g', norm( ErrEv, 1 ), norm( ErrEv, 2 ), norm( ErrEv, 'inf' ) ) );
    disp( sprintf( 'rel norms: 1norm: %12g,  2norm: %12g,  infnorm: %12g', norm(ErrEv,1)/norm(Ev,1), norm(ErrEv,2)/norm(Ev,2), norm(ErrEv,'inf')/norm(Ev,'inf') ) );
end;