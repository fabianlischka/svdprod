% $Id$

wtf = 0;
if wtf
    fid = fopen( 'symtridhh.log', 'w' );
else
    fid = 1;
end;

tol = 1e-10;
fprintf( fid,  'TEST SYMTRIDHH\n' );
fprintf( fid,  'We generate symmetrized testmatrices A, and compute tridiag T = Q''AQ,\nwith Q orthogonal. ' );
fprintf( fid,  'We then compute evals of A and T, and display\nthe norm of the difference (relative to the' );
fprintf( fid,  'norm of the evals itself).\nWe also consider the norm of redidual A-QTQ'', and check\n' );
fprintf( fid,  'orthogonality, namely the norm of Q''Q-I.\nAll errors should be around 1e-15.\n' );
for N=3:37:77
    fprintf( fid,  '\nTesting SYMTRIDHH with various (symmetrized) %g x %g matrices\n', N, N );
    fprintf( fid,  'Type:   Eig 2-norm        RelRes          Orth\n' );
    for Typ = 1:14
        TestMat = gentestmat( Typ, N, N, 0 );
        TestMat = TestMat + TestMat';
        Ev      = eig( TestMat );
        [ T, Q ]= symtridhh( TestMat );
        ErrEv   = norm( sort( eig(T) ) - sort( Ev ) ) / norm( Ev );
        Res     = norm( TestMat - Q*T*Q' );
        RelRes  = Res / norm( TestMat );
        Orth    = norm( Q'*Q - eye( N ) );
        fprintf( fid,  '  %2g: %12g, %12g, %12g\n', Typ, ErrEv, RelRes, Orth );
        if RelRes > tol || Orth > tol || ErrEv > tol
            fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
        end
    end;
end;

if wtf; fclose( fid ); end;
