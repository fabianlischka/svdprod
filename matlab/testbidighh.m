% $Id$

N = 30;
for M = 30:21:72
    disp( sprintf( '\nTesting BIDIGHH with various %g x %g matrices', M, N )); 
    for Typ = 1:14
        TestMat = gentestmat( Typ, M, N, 0 );
        for trans = 0:(M==N)
            [B,U,V] = bidighh( TestMat );
            % B       = [ diag( BS(:,1) ) + diag( BS(1:end-1,2), 1 ); zeros( M-N, N ) ];
            dev     = norm( B - U'*TestMat*V );
            reldev  = dev / norm( TestMat );
            svdB    = sort( svd( B ) );
            svdA    = sort( svd( TestMat ) );
            devsvd  = norm( svdA - svdB, 'inf' );
            Udev    = norm( U'*U - eye( M ) );
            Vdev    = norm( V'*V - eye( N ) );
            disp( sprintf( 'type %2g: err: %12g, relerr: %12g, relerrSVs: %12g, Udev: %12g, Vdev: %12g', Typ, dev, reldev, devsvd/norm( svdA, 'inf' ), Udev, Vdev ) );
            TestMat = TestMat';
        end;
    end;
end;
