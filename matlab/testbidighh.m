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
    for Typ = 1:14
        TestMat = gentestmat( Typ, M, N, 0 );
        for trans = 0:(M==N)
            [BS,U,V]= bidighh( TestMat );
            B       = [ diag( BS(:,1) ) + diag( BS(1:end-1,2), 1 ); zeros( M-N, N ) ];
            dev     = norm( B - U'*TestMat*V );
            reldev  = dev / norm( TestMat );
            svdB    = sort( svd( B ) );
            svdA    = sort( svd( TestMat ) );
            devsvd  = norm( svdA - svdB, 'inf' )/norm( svdA, 'inf' );
            Udev    = norm( U'*U - eye( M ) );
            Vdev    = norm( V'*V - eye( N ) );
            fprintf( fid, 'type %2g: err: %12g, relerr: %12g, relerrSVs: %12g, Udev: %12g, Vdev: %12g\n', Typ, dev, reldev, devsvd, Udev, Vdev );
            if reldev > tol || Udev > tol || Vdev > tol || devsvd > tol
                fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
            end
            TestMat = TestMat';
        end;
    end;
end;

if wtf; fclose( fid ); end;
