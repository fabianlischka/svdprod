% $Id$

wtf = 0;
if wtf
    fid = fopen( 'symtridhh.log', 'w' );
else
    fid = 1;
end;

tol = 1e-10;
fprintf( fid,  '\n\nTEST SYMTRIDHH\n' );
fprintf( fid,  'We generate symmetrized testmatrices A, and compute tridiag T = Q''AQ, with Q orthogonal.\n' );
fprintf( fid,  'We then compute evals of A and T, and display the norm of the difference (relative to the\n' );
fprintf( fid,  'norm of the evals itself). We also consider the norm of redidual A-QTQ'', and check\n' );
fprintf( fid,  'orthogonality, norm of Q''Q-I. All errors should be around 1e-15.\n' );
for N=3:37:77
    fprintf( fid,  '\nTesting SYMTRIDHH with various (symmetrized) %g x %g matrices\n', N, N ); 
    for Typ = 1:14
        TestMat = gentestmat( Typ, N, N, 0 );
        TestMat = TestMat + TestMat';
        Ev      = eig( TestMat );
        [ T, Q ]= symtridhh( TestMat );
        ErrEv   = norm( sort( eig(T) ) - sort( Ev ) ) / norm( Ev );
        Res     = norm( TestMat - Q*T*Q' );
        RelRes  = Res / norm( TestMat );
        Orth    = norm( Q'*Q - eye( N ) );
        fprintf( fid,  'Type %2g: eval infnorm: %12g; RelRes: %12g, Orth: %12g\n', Typ, ErrEv, RelRes, Orth ) ;
        if RelRes > tol || Orth > tol || ErrEv > tol
            fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
        end
    end;
end;

if wtf; fclose( fid ); end;
