
for L = 1:3
    Liter   = 3;
Typ     = err(1,1,L);
ts      = sprintf( 'Relative Error: Type %g, ', Typ );

    K           = err(2,1,L);
    N           = err(3,1,L);
    implerr     = abs( err(1:N,3,L)-err(1:N,2,L) );
    explerr     = abs( err(1:N,4,L)-err(1:N,2,L) );
    true1       = err(1:N,2,L);
    implrel1    = err(1:N,5,L);
    explrel1    = err(1:N,6,L);
    ts          = [ ts sprintf( '%g matrices %gx%g (-)\n', K, N, N ) ];
    L = L + Liter;
    K           = err(2,1,L);
    N           = err(3,1,L);
    implerr     = abs( err(1:N,3,L)-err(1:N,2,L) );
    explerr     = abs( err(1:N,4,L)-err(1:N,2,L) );
    true2       = err(1:N,2,L);
    implrel2    = err(1:N,5,L);
    explrel2    = err(1:N,6,L);
    ts          = [ ts sprintf( '%g matrices %gx%g (--), ', K, N, N ) ];    
    L = L + Liter;
    K           = err(2,1,L);
    N           = err(3,1,L);
    implerr     = abs( err(1:N,3,L)-err(1:N,2,L) );
    explerr     = abs( err(1:N,4,L)-err(1:N,2,L) );
    true3       = err(1:N,2,L);
    implrel3    = err(1:N,5,L);
    explrel3    = err(1:N,6,L);
    ts          = [ ts sprintf( '%g matrices %gx%g, (:)', K, N, N ) ];
    clf reset;
     loglog( true1, implrel1, '-+', true1, explrel1, '-o', ...
            true2, implrel2,'--+', true2, explrel2,'--o', ...
            true3, implrel3, ':+', true3, explrel3, ':o'  ...
        );
    lol = gcf;
%    semilogy( 1:N, implrel, '-+', 1:N, explrel, '-o' );
    title(ts);
    if Liter == 1
        name = sprintf( 're%1gK%gNNN', Typ, K );
    else
        name = sprintf( 're%1gKKKN%g', Typ, N );
    end;
    name
    saveas( lol, name, 'eps' )
end;
end;