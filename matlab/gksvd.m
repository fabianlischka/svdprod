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
        if B( k, 1 ) == 0   % FIX FIX: compare with some eps???
            err;            % rotate that row away, and look again
        else    % do some real work
            [ B, U, V ] = gksvdstep( B, U, V );
        end;    
    end;
end;