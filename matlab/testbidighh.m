% $Id$

N = 20;
for M = 20:10:30
    disp( sprintf( '\nTesting BIDIGHH with various %g x %g matrices', M, N )); 
    for Typ = 1:14
        TestMat = gentestmat( Typ, M, N, 0 );
        disp( sprintf( 'for type %g', Typ ) );
        for trans = 0:(M==N)
            [B,U,V] = bidighh( TestMat );
            dev     = norm( B - U'*TestMat*V );
            reldev  = dev / norm( TestMat );
            svdB    = sort( svd( B ) );
            svdA    = sort( svd( TestMat ) );
            devsvd  = norm( svdA - svdB, 'inf' );
            disp( sprintf( 'err: %12g,  rel err: %12g,  norm SVs: %12g, rel norm SVs: %12g', dev, reldev, devsvd, devsvd/norm( svdA, 'inf' ) ) );
            TestMat = TestMat';
        end;
    end;
end;
