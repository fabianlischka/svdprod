% function [ D, U, V ] = gksvd( A, tol )
% GKSVD computes the SVD of A
% Given A MxM and tol, this computes U, V MxM orthogonal and D diagonal 
% such that A=U(D+E)V', where two norm E is around unit roundoff x two norm A
% note: will later extend this to M >= N

% BIDIGHH computes bidiagonal B=U'AV using householder reflections
[ B, U, V ] = bidighh( A );

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
if M < N
    error( 'A must have M >= N' );
end;

q = 0;
while q < N
    q   = 0; p = N-1;
    for k = N-1:-1:1
        if abs( B(k,2) ) < tol*( abs( B(k,1) ) + abs( B(k+1,1) ) )
            B(k,2) = 0;
            if q == N-k-1   % state one
                q = N-k;
                p = k-1;
            % else          % state three
            end;
        else
            if p == k       % state two
                p = k-1;
            % else          % state three
            end;
        end;
    end;
    if q == N-1
        q=N;    % and abort...
    else        % do some work
        k = p+1;
        while B( k, 1 ) ~= 0  && k < N-q
            k = k+1;
        end;
        % now, k == N-q, or B(k,1) == 0, or both
        if abs( B( k, 1 ) ) <= eps * norm( B, 'inf' )
            % now, B( k, 1 ) == 0, p+1 <= k <= N-q, B(N-q,2) == 0 (note:
            % even if q == 0, by the format we have chosen)
            % So, to do: if k<N-q, rotate that row away with Givens rotations from
            % left, G(k,j), j=k+1:N-q
            % if k==N-q, then apply Givens from the 
            % right, G(k,j), j=N-q-1:-1:p+1
            disp('rotating empty diagonal away')
            if k < N-q
                bulge = B(k,2);
                B(k,2)= 0;
                for j=k+1:N-q
                    [ c s ] = givens( B(j,1), bulge );
                    % need to transpose (since mult from left), and 
                    % transpose again, since need to eliminate upper
                    % element
                    G     = [ c s; -s c ]; 
                    % rotate B from left: B+ = G B,  G=G(k,j)
                    % zerod =  c * bulge + s * B(j,1)
                    B(j,1)= -s * bulge + c * B(j,1);
                    bulge =  s * B(j,2);
                    B(j,2)=  c * B(j,2);                    
                    % rotate U from right (since transposed): U+ = U G'
                    U(:,[k j]) = U(:,[k j]) * G';       % flops: 6*M
                end
            else
                for j=N-q-1:-1:p+1
                    err;
                    [ c s ] = givens( 1,2 );
                    rotate B from right
                    rotate V from right
                end
            end
        else    % do some real work - GK SVD step for p+1:N-q
            [ B(p+1:N-q,:), U(:,p+1:N-q), V(:,p+1:N-q) ] = gksvdstep( B(p+1:N-q,:), U(:,p+1:N-q), V(:,p+1:N-q) );
        end;    
    end;
    % B(:,2)'
    % disp( sprintf( '%12g   %12g   %12g', (B(1,2)), log(B(2,2)), log(B(3,2)) ) );
end;
D = abs( B( :,1 ) );