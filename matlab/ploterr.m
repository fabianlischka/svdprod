
for L=1:  size(err,3)
    Typ         = err(1,1,L);
    K           = err(2,1,L);
    N           = err(3,1,L);
    implerr     = abs( err(1:N,3,L)-err(1:N,2,L) );
    explerr     = abs( err(1:N,4,L)-err(1:N,2,L) );
    true        = err(1:N,2,L);
    implrel     = err(1:N,5,L);
    explrel     = err(1:N,6,L);
%    loglog( true, implrel, '-+', true, explrel, '-o' );
    semilogy( 1:N, implrel, '-+', 1:N, explrel, '-o' );
    title(sprintf('Relative error: Type %g, prod of %g matrices, %g x %g',Typ,K,N,N));
    pause
end;