% $Id$

tol = 1e-10;

disp( sprintf( '\n\nTEST GKSVD' ) );
disp( 'We compute singular values of a test matrix using svd(), and then using our implementation of' );
disp( 'the Golub-Kahan SVD alg, and display the norm of the difference, relative to the norm of the svs.' );
disp( 'We also test the norm of A-UDV'', and the orthogonality, via norm U''U-I and V''V-I.' );
disp( 'All errors should be around 1e-15.' );
for M=5:21:68
    N=min(30,M);
    disp( sprintf( '\nTesting GKSVD with various %g x %g matrices', M, N )); 
    for Typ = 1:14
        TestMat = gentestmat( Typ, M, N, 0 );
        [D,U,V] = gksvd( TestMat, 1e-14 );
        svds    = svd( TestMat );
        ErrSvd  = norm( sort( D ) - sort( svds ) ) / norm( svds );
        Res     = norm( TestMat - U(:,1:N)*diag(D)*V' );
        RelRes  = Res / norm( TestMat );
        UOrth   = norm( U'*U - eye( M ) );
        VOrth   = norm( V'*V - eye( N ) );        
        disp( sprintf( 'Type %2g: svd 2-norm: %12g, RelRes: %12g, UOrth: %12g, VOrth: %12g', Typ, ErrSvd, RelRes, UOrth, VOrth ) );
        if RelRes > tol || UOrth > tol || VOrth > tol || ErrSvd > tol
            disp( '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^' );
        end
    end;
end

N = 4;
disp( sprintf( '\nTesting GKSVD, all combinations of off-diagonal elements zero, %g x %g matrices', N, N )); 
for Typ=0:(2^N-1)
    Bits    = bitget( Typ, 1:N-1 );
    TestMat = -diag( 1:N ) + diag( Bits, 1 );
        [D,U,V] = gksvd( TestMat, 1e-14 );
        svds    = svd( TestMat );
        ErrSvd  = norm( sort( D ) - sort( svds ) ) / norm( svds );
        Res     = norm( TestMat - U*diag(D)*V' );
        RelRes  = Res / norm( TestMat );
        UOrth   = norm( U'*U - eye( N ) );
        VOrth   = norm( V'*V - eye( N ) );        
        disp( sprintf( 'Type %2g: svd 2-norm: %12g, RelRes: %12g, UOrth: %12g, VOrth: %12g', Typ, ErrSvd, RelRes, UOrth, VOrth ) );
        if RelRes > tol || UOrth > tol || VOrth > tol || ErrSvd > tol
            disp( '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^' );
            err;
        end
end;
