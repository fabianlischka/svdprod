function A = cholesky( A );
% CHOLESKY finds the cholesky decomposition of A=F'F, 
% F uppertriangular, and returns F

% reference: Golub, Van Loan; 3rd ed; ch. 4.2.5
% $Id$
N = size(A,1);

% outer product version
for k=1:N
    if A(k,k) < 0
        error( 'A not positive definite' );
    end;
    A(k,k)      = sqrt( A(k,k) );
    A(k,k+1:end)= A(k,k+1:end) / A(k,k);
                        % flops per iter: 2(N-k)
                        % total flops: N^2
    for j=k+1:N
        A(j,j:end) = A(j,j:end) - A(k,j)*A(k,j:end);
                        % flops per iter: 2(N-j)
                        % flops per k: (N-k)^2
                        % total flops: N^3/3
    end;
    A(k+1:N,k) = 0;
end;
% total flop count: 1/3 N^3
