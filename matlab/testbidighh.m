% $Id$

N = 20;
for M = 20:10:30
    disp( sprintf( '\nTesting BIDIGHH with various %g x %g matrices', M, N )); 
    for Typ = 1:13
        TestMat = gentestmat( Typ, N, N, 0 );
        disp( sprintf( 'for type %g', Typ ) );
        for trans = 0:1
            [B,U,V] = bidighh( TestMat );
            dev     = norm( B - U'*TestMat*V );
            reldev  = dev / norm( TestMat );
            svdB    = sort( svd( B ) );
            svdA    = sort( svd( TestMat ) );
            disp( sprintf( 'err norm: %12g,  rel norm: %12g,  infnorm SVs: %12g', dev, reldev, norm( svdA - svdB, 'inf' ) ) );
            TestMat = TestMat';
        end;
    end;
end;
