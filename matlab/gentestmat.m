function A = gentestmat( k, M, N, verbose )
% GENTESTMAT( k, M, N ) generates a MxN testmatrix of type k (M>=N)
% $Id$

As      = { '-triu( ones( M ) ) + 2 * diag( ones( M, 1 ) )', 
                'uppertriang, structured, eval = 1, only 1 evect, bad condition';
            'triu( 10*rand( M ) - 1 ) - 20 * diag( ones( M, 1 ) )',
                'uppertriang, large negative diag';
            '2*rand( M, N ) - 1',
                'elements uniformly random -1..1';
            'exp(16*rand( M, N )-15).*(rand( M, N )-.5)',  
                'elements of very different size';
            'U* diag(2.^(-1:-1:-N)) *V',         
                'very different singular values: 1/2 .. 1/(2^N)';
            '10 * randn( M, 1 ) * rand( 1, N ) + 1e-15*rand( M, N )', 
                'rank 1 matrix with slight perturbations';
            'W* ( diag(1e-6*ones(M,1)) + diag(ones(M-1,1),1) ) * transpose(W)',
                '1 small eval, one evect';
            'W* ( diag(1e+6*ones(M,1)) + diag(ones(M-1,1),1) ) * transpose(W)',
                '1 big   eval, one evect';
            'W* diag([1:M-1], 1 ) * transpose(W)',
                'zero eigenvals, one evect',
            '-(rand( M, N )>(1-5/(M*N)))',
                'sparse matrix with 95% negative ones, rest zero'
            'hilb( M )',
                'Hilbert matrix';
            'full( gallery( ''dorr'', M, 0.02 ) )',
                'Dorr matrix';
            'gallery(''kahan'', [ M, N ] )'
                'Kahan matrix';
          };

if nargin == 2 || N > M
    N = M; 
end;
if nargin < 4
    verbose = 1;
end;

% produce some random orthogonal matrices, needed for some types
[U,X]=qr( rand( M, N ), 0 );
[V,X]=qr( rand( N ) );
[W,X]=qr( rand( M ) );

if verbose
disp( sprintf( '\nGenerating %d x %d matrix, type %d (%s): generated by\n %s:\n', M, N, k, As{k*2}, As{ k*2-1 } ) );
end;
A = eval( As{ k*2-1 } );
A = A(1:M, 1:N);
