function [ D, Q ] = symqrschur( A, tol )
% SYMQRSCHUR computes an approximate sym Schur decomposition Q'AQ = D
% given symmetric real A, and a tolerance tol,
% using implied QR steps with Wilkinson shift

% reference: Golub, Van Loan; 3rd ed; ch. 8.3.5, alg. 8.3.3
% $Id$

% T = Q'*A*Q is (sym) tridiagonal
[ T, Q ] = symtridhh( A );
if nargin == 1
    tol = 1e-14;
end;

TS  = stridl2s( T );  % TS(k,1) == T(k,k); TS(k,2) == T(k,k+1) == T(k+1,k)
N   = size( TS, 1 );
q   = 0;
while q < N
    q   = 0; p = N-1;
    % disp(' - ');
    for k = N-1:-1:1
        % disp( sprintf( 'off diag: %12g, on diags: %12g   %12g', abs( TS( k, 2 ) ), abs( TS( k, 1 ) ), abs( TS( k+1, 1 ) ) ) );
        if abs( TS( k, 2 ) ) <= tol * ( abs( TS( k,1 ) ) + abs( TS( k+1, 1 ) ) )
            TS( k, 2 ) = 0;
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
    % after this loop: TS(p+1,2)..TS(N-q-1,2) are non zero,
    % TS(N-q,2)..TS(N-1,2) are zero
    % ie all zero -> q=N, p=0; all nonzero -> q=0; p=0
    % three states: 3, 2, 1
    
    % start: q = 0; p = N-1; k = N-1
    % state 1:
    % TS(k,2) == 0 ? q++, k--, p--, state 1  :  k--, state 2
    % state 2:
    % TS(k,2) == 0 ? k--, state 3 : k--, p--, state 2
    % state 3:
    % k--
    
    if q == N-1
        q = N;  % and abort
    else        % do some work: on D22, ie TS( p+1:N-q, : )
        [ TS( p+1:N-q, : ), Q( p+1:N-q, p+1:N-q ) ] = qrimstep( TS( p+1:N-q, : ), Q( p+1:N-q, p+1:N-q ) );
    end;
%    TS
%    Q
%    Q*strids2l(TS)*Q'
end;
D = TS( :,1 );