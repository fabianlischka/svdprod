% $Id$

N = 60;

disp( sprintf( '\nTesting SYMQRSCHUR with various (symmetrized) %g x %g matrices', N, N )); 
for Typ = 1:12
	TestMat = gentestmat( Typ, M, N, 0 );
    TestMat = TestMat + TestMat';
    Ev      = eig( TestMat );
    ErrEv   = sort(symqrschur( TestMat ) ) - sort( Ev );
    disp( sprintf( 'for type %g', Typ ) );
    disp( sprintf( 'err norms: 1norm: %12g,  2norm: %12g,  infnorm: %12g', norm( ErrEv, 1 ), norm( ErrEv, 2 ), norm( ErrEv, 'inf' ) ) );
    disp( sprintf( 'rel norms: 1norm: %12g,  2norm: %12g,  infnorm: %12g', norm(ErrEv,1)/norm(Ev,1), norm(ErrEv,2)/norm(Ev,2), norm(ErrEv,'inf')/norm(Ev,'inf') ) );
end;