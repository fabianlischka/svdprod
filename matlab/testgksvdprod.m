function testgksvdprod
% $Id$

global err;
wtf = 0;
if wtf
    fid = fopen( 'gksvdprod.log', 'w' );
else
    fid = 1;
end;
tol = 1e-10;

fprintf( fid,  'TEST GKSVDPROD\n' );

if 0 % generate detailed error data
Layer = 1;
for Typ = 1:2
    for K = 5:8:21
        for N = 8:12:32
            fprintf( fid,  'computing Type %g, %3g matrices %4g x %4g\n', Typ, K, N, N );
            [ A, true ] = gentestprod( Typ, N, K );
            impl    = sort( gksvdprod(A, 1e-16) );
			Amul    = eye(size(A,1)); for k=1:size(A,3); Amul = A(:,:,k)*Amul; end; 
			expl    = sort( svd(Amul) );
            implrel = abs( (impl - true) ./ true );
            explrel = abs( (expl - true) ./ true );
            % err Nx5 contains in col 1 [Typ K N 0000000]'
            % col 2 true, col 3 impl, 4 expl, col 5 implrel, 6 explrel
            err(1:3,1,  Layer)  = [Typ;K;N];
            err(1:N,2:6,Layer)  = [true impl expl implrel explrel];
            Layer = Layer + 1;
        end
    end
end
end

K = 4;
for M=5:21:47 % 47
    fprintf( fid,  '\nTesting GKSVDPROD with various %g x %g matrices\n', M, M ); 
    fprintf( fid,  'Prod of types:   SVs 2-norm        RelRes        U Orth        V Orth\n' );
    TestMats = zeros(M,M,14);
    for Typ = 1:14
        TestMats(:,:,Typ) = gentestmat( Typ, M, M, 0 );
    end;
    for k=1:(14-K)
        A       = TestMats(:,:,k:k+K);
        [D,U,V] = gksvdprod( A, 1e-16 );
        Amul    = eye(size(A,1)); for kk=1:size(A,3); Amul = A(:,:,kk)*Amul; end; 
        expl    = sort( svd( Amul ) );
        impl    = sort( D );
        ErrSvd  = norm( impl - expl ) / norm( expl );
        Res     = norm( Amul - U*diag(D)*V' );
        RelRes  = Res / norm( Amul );
        UOrth   = norm( U'*U - eye( M ) );
        VOrth   = norm( V'*V - eye( M ) );        
        fprintf( fid,  '   %3g ... %2g: %12g, %12g, %12g, %12g\n', k, k+K, ErrSvd, RelRes, UOrth, VOrth );
        if RelRes > tol || UOrth > tol || VOrth > tol || ErrSvd > tol
            fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
        end
    end;
end

M = 8;
K = 6;
A   = zeros( M,M,K );
fprintf( fid,  '\nTesting GKSVDPROD, random combinations of off-diagonal elements zero, %g x %g matrices\n', M, M ); 
fprintf( fid,  'Combinations:   SVs 2-norm        RelRes        U Orth        V Orth\n' );
for l=1:10
    for k = 1:K
        Bits    = bitget( floor( rand(1) * 2^M ), 1:M-1 );
        A(:,:,k)= -diag( 1:M ) + diag( Bits, 1 );
    end;
    [D,U,V] = gksvdprod( A, 1e-16 );
    Amul    = eye(size(A,1)); for kk=1:size(A,3); Amul = A(:,:,kk)*Amul; end; 
    expl    = sort( svd( Amul ) );
    impl    = sort( D );
    ErrSvd  = norm( impl - expl ) / norm( expl );
    Res     = norm( Amul - U*diag(D)*V' );
    RelRes  = Res / norm( Amul );
    UOrth   = norm( U'*U - eye( M ) );
    VOrth   = norm( V'*V - eye( M ) );        
    fprintf( fid,  '       %2g   : %12g, %12g, %12g, %12g\n', K, ErrSvd, RelRes, UOrth, VOrth );
    if RelRes > tol || UOrth > tol || VOrth > tol || ErrSvd > tol
        fprintf( fid,  '^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^\n' );
    end
end;

function [ A, true ] = gentestprod( Typ, N, K )
A = zeros(N,N,K);
switch Typ
    case 1
        S = diag(2.^(-1:-1:-N));
        [V,X]=qr(randn(N));
        % AAAclose = V;
		for k=1:K
            [U,X]=qr(randn(N));
            A(:,:,k) = U'*S*V;
            V=U;
		end;
        true = 2.^(K*(-1:-1:-N))';
        % AAAclose = U'*diag(true)*V;
    case 2
        a = 2;
        b = -1;
        Toep = toeplitz( [ a b zeros(1,N-2) ] );
        for k=1:K
            A(:,:,k) = Toep;
		end;
		true = (a + 2*cos((1:N)'/(N+1) * pi )).^K;
end;
true = sort( true );

if wtf; fclose( fid ); end;
