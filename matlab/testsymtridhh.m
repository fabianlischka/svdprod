% $Id$

tol = 1e-10;
disp( sprintf( '\n\nTEST SYMTRIDHH' ) );
disp( 'We generate symmetrized testmatrices A, and compute tridiagonal T = Q''AQ, with Q orthogonal.' );
disp( 'We then compute evals of A and T, and display the norm of the difference (relative to the' );
disp( 'norm of the evals itself). We also consider the norm of redidual A-QTQ'', and check' );
disp( 'orthogonality, norm of Q''Q-I. All errors should be around 1e-15.' );
for N=3:37:77
    disp( sprintf( '\nTesting SYMTRIDHH with various (symmetrized) %g x %g matrices', N, N )); 
    for Typ = 1:14
        TestMat = gentestmat( Typ, N, N, 0 );
        TestMat = TestMat + TestMat';
        Ev      = eig( TestMat );
        [ T, Q ]= symtridhh( TestMat );
        ErrEv   = norm( sort( eig(T) ) - sort( Ev ) ) / norm( Ev );
        Res     = norm( TestMat - Q*T*Q' );
        RelRes  = Res / norm( TestMat );
        Orth    = norm( Q'*Q - eye( N ) );
        disp( sprintf( 'Type %2g: eval infnorm: %12g; RelRes: %12g, Orth: %12g', Typ, ErrEv, RelRes, Orth ) );
        if RelRes > tol || Orth > tol || ErrEv > tol
            disp( '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^' );
        end
    end;
end;