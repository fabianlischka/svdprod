% $Id$

wtf = 0;
if wtf
    fid = fopen( 'symqrschur.log', 'w' );
else
    fid = 1;
end;
tol = 1e-10;

fprintf( fid,  '\n\nTEST SYMQRSCHUR\n' );
fprintf( fid,  'We compute eigenvalues of a test matrix using eig(), and then using the schur decomposition\n' );
fprintf( fid,  'computed by symqrschur, and display the norm of the difference of the vector of eigenvalues\n' );
fprintf( fid,  'relative to the norm of the vector of eigenvalues itself. Also the norm of A-QDQ'', and Q''Q-I.\n' );
fprintf( fid,  'All errors should be around 1e-15.\n' );
for N=5:25:55
    fprintf( fid,  '\nTesting SYMQRSCHUR with various (symmetrized) %g x %g matrices\n', N, N ); 
    for Typ = 1:14
        TestMat = gentestmat( Typ, N, N, 0 );
        TestMat = TestMat + TestMat';
        Ev      = eig( TestMat );
        [ D,Q ] = symqrschur( TestMat );
        ErrEv   = norm( sort( D ) - sort( Ev ) )/norm( Ev );
        Res     = norm( TestMat - Q*diag(D)*Q' );
        RelRes  = Res / norm( TestMat );
        Orth    = norm( Q'*Q - eye( N ) );
        fprintf( fid,  'Type %2g: eval 2-norm: %12g, RelRes: %12g, Orth: %12g\n', Typ, ErrEv, RelRes, Orth );
        if RelRes > tol || Orth > tol || ErrEv > tol
            fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
        end
    end;
end

N = 4;
fprintf( fid,  '\nTesting SYMQRSCHUR, all combinations of off-diagonal elements zero, %g x %g matrices\n', N, N ); 
for Typ=0:(2^N-2)
    Bits    = bitget( Typ, 1:N-1 );
    TS      = [ 1:N; Bits, 0 ]';
    TestMat = strids2l( TS );
    Ev      = eig( TestMat );
    [ D,Q ] = symqrschur( TestMat );
    ErrEv   = norm( sort( D ) - sort( Ev ) )/norm( Ev );
    Res     = norm( TestMat - Q*diag(D)*Q' );
    RelRes  = Res / norm( TestMat );
    Orth    = norm( Q'*Q - eye( N ) );
    fprintf( fid,  'Type %2g: eval 2-norm: %12g, RelRes: %12g, Orth: %12g\n', Typ, ErrEv, RelRes, Orth );
    if RelRes > tol || Orth > tol || ErrEv > tol
        fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
    end
end;

if wtf; fclose( fid ); end;
