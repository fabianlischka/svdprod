% $Id$

wtf = 0;
if wtf
    fid = fopen( 'gksvd.log', 'w' );
else
    fid = 1;
end;
tol = 1e-10;

fprintf( fid,  '\n\nTEST GKSVD\n' );
fprintf( fid,  'We compute singular values of a test matrix using svd(), and then using our implementation of\n' );
fprintf( fid,  'the Golub-Kahan SVD alg, and display the norm of the difference, relative to the norm of the svs.\n' );
fprintf( fid,  'We also test the norm of A-UDV'', and the orthogonality, via norm U''U-I and V''V-I.\n' );
fprintf( fid,  'All errors should be around 1e-15.\n' );
for M=5:21:68
    N=min(30,M);
    fprintf( fid,  '\nTesting GKSVD with various %g x %g matrices\n', M, N );
    fprintf( fid,  'Type:   svd 2-norm        RelRes        U Orth        V Orth\n' );
    for Typ = 1:14
        TestMat = gentestmat( Typ, M, N, 0 );
        [D,U,V] = gksvd( TestMat, 1e-14 );
        svds    = svd( TestMat );
        ErrSvd  = norm( sort( D ) - sort( svds ) ) / norm( svds );
        Res     = norm( TestMat - U(:,1:N)*diag(D)*V' );
        RelRes  = Res / norm( TestMat );
        UOrth   = norm( U'*U - eye( M ) );
        VOrth   = norm( V'*V - eye( N ) );        
        fprintf( fid,  '  %2g: %12g, %12g, %12g, %12g\n', Typ, ErrSvd, RelRes, UOrth, VOrth );
        if RelRes > tol || UOrth > tol || VOrth > tol || ErrSvd > tol
            fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
        end
    end;
end

N = 4;
fprintf( fid,  '\nTesting GKSVD, all combinations of off-diagonal elements zero, %g x %g matrices\n', N, N );
fprintf( fid,  'Type:   svd 2-norm        RelRes        U Orth        V Orth\n' );
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
        fprintf( fid,  '  %2g: %12g, %12g, %12g, %12g\n', Typ, ErrSvd, RelRes, UOrth, VOrth );
        if RelRes > tol || UOrth > tol || VOrth > tol || ErrSvd > tol
            fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
            err;
        end
end;

if wtf; fclose( fid ); end;
