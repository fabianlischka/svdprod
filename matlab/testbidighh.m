% $Id$

wtf = 0;
if wtf
    fid = fopen( 'bidighh.log', 'w' );
else
    fid = 1;
end;

tol = 1e-12;
N   = 30;
for M = 30:21:72
    fprintf( fid, '\nTesting BIDIGHH with various %g x %g matrices\n', M, N );
    fprintf( fid, 'Type:  B-U''AV Res  Relative Res   SVs Rel Err        U Orth        V Orth\n' );
    for Typ = 1:14
        TestMat = gentestmat( Typ, M, N, 0 );
        [BS,U,V]= bidighh( TestMat );
        B       = [ diag( BS(:,1) ) + diag( BS(1:end-1,2), 1 ); zeros( M-N, N ) ];
        dev     = norm( B - U'*TestMat*V );
        reldev  = dev / norm( TestMat );
        svdB    = sort( svd( B ) );
        svdA    = sort( svd( TestMat ) );
        devsvd  = norm( svdA - svdB, 'inf' )/norm( svdA, 'inf' );
        Udev    = norm( U'*U - eye( M ) );
        Vdev    = norm( V'*V - eye( N ) );
        fprintf( fid, ' %2g: %12g, %12g, %12g, %12g, %12g\n', Typ, dev, reldev, devsvd, Udev, Vdev );
        if reldev > tol || Udev > tol || Vdev > tol || devsvd > tol
            fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
        end
    end;
end;

if wtf; fclose( fid ); end;
