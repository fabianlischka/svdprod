function [ D, U, V ] = gksvdsteps( B, U, V, tol )
% GKSVDSTEPS computes the SVD of a bidiagonal B, and
% overwrites U, V such that if previously B=U'AV, now D=U'AV
% this is an internal function, so no errorchecking on input
% both B and D in short format (as column vectors)

% reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.2
% $Id$

M = size( U, 1 );   % rows
N = size( V, 1 );   % cols

% partition B into an unreduced part in the middle, and a diagonal part at
% the bottom
q = 0;
while q < N
    q   = 0; p = N-1;
    for k = N-1:-1:1
        % note: below, we need <=, not <, otherwise we get caught in an
        % infitite loop, if these elements are exactly zero!
        if abs( B(k,2) ) <= tol*( abs( B(k,1) ) + abs( B(k+1,1) ) )
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
        % NOW, with D11 = BL(1:p,1:p), D22 = BL(p+1:N-q,p+1:N-q),
        % D33 = BL(N-q+1:N, N-q+1:N), D33 is diagonal, D22 unreduced
    else        % do some work
        % first, determine if any element on the diagonal is zero
        % if so, rotate it aside
        k = p+1;
        smalldiag = tol * norm( B, 'inf' );
        while abs( B( k, 1 ) ) > smalldiag  && k < N-q
            k = k+1;
        end;
        % now, k == N-q, or abs( B(k,1) ) <= smalldiag, or both
        if abs( B( k, 1 ) ) <= smalldiag
            B( k, 1 ) = 0;
            % now, B( k, 1 ) approx 0, p+1 <= k <= N-q, B(N-q,2) == 0 (note:
            % even if q == 0, by the format we have chosen)
            if k < N-q
                % to do: if k<N-q, rotate that row away with Givens rotations from
                % left, G(k,j), j=k+1:N-q
                bulge = B(k,2);
                B(k,2)= 0;
                for j=k+1:N-q
                    [ c s ] = givens( B(j,1), bulge );
                    % need to transpose (since mult from left), and 
                    % transpose again, since need to eliminate upper
                    % element
                    % rotate B from left: B+ = G B,  G=G(k,j)
                    % zerod =  c * bulge + s * B(j,1)
                    B(j,1)  = -s * bulge + c * B(j,1);
                    bulge   =  s * B(j,2);
                    B(j,2)  =  c * B(j,2);                    
                    % rotate U from right (since transposed): U+ = U G'
                    G       = [ c s; -s c ]; 
                    U(:,[k j]) = U(:,[k j]) * G';       % flops: 6*M
                end % for j
            else % k==N-q
                % if k==N-q, then apply Givens from the 
                % right, G(k,j), j=N-q-1:-1:p+1
                bulge = B(N-q-1,2);
                B(N-q-1,2)= 0;
                for j=N-q-1:-1:p+1
                    % rotate B from right, B+ = B G;  G=G(k,j)
                    [ c s]  = givens( B(j,1), bulge );
                    % zerod   = s * B(j,1) + c * bulge
                    B(j,1)  = c * B(j,1) - s * bulge;
                    if j>p+1
                        bulge   = s * B(j-1,2);
                        B(j-1,2)= c * B(j-1,2);
                    end;
                    % rotate V from right
                    G       = [ c s; -s c ]; 
                    V(:,[j k]) = V(:,[j k]) * G;       % flops: 6*M
                end
            end % if k
        else    % no element on the diagonal is zero
            % do some real work - GK SVD step for p+1:N-q
            [ B(p+1:N-q,:), U(:,p+1:N-q), V(:,p+1:N-q) ] = gksvdstep( B(p+1:N-q,:), U(:,p+1:N-q), V(:,p+1:N-q) );
        end; 
    end; % if q
end; % while q<N

% now fix sign - want cols of V such that elements of D are non-negative
D = B( :,1 );
for k = 1:N
    if D(k) < 0
        D(k)    = -D(k);
        V(:,k)  = -V(:,k);
    end;
end;