function [ svals, t2series ] = svals( At, Params )
% SVALS computes the singular values of At by applying two sided Jacobi (=Givens)
% rotations to reduce the size of off-diagonal elements after symmetrizing
% then returns the absolute value of the diagonal

% reference: Golub, Van Loan; 3rd ed; ch. 8.4, ch. 8.6.3
% $Id$


if ~isfield( Params,'fun');         Params.fun      = 'maxsumsq';       end;
if ~isfield( Params,'maxiter');     Params.maxiter  = 100;              end;
if ~isfield( Params,'dispiter');    Params.dispiter = 20;               end;

% for the time being, we deal with non-square matrices by obtaining a QR
% factorization, and working only with the uppertriangualr R (which has
% zeros below and is effectively square)
if size( At, 1 ) ~= size( At, 2 )
    [ Qt, At ] = qr( At, 0 );
end;

% we terminate when tsquare does not improve anymore
t2       = norm( At, 'fro') - norm( diag(At),'fro');   
t2series = [ t2 * 10 ];         % initial value for first iteration
iter     = 0;

while t2 < t2series( end:end ) && iter < Params.maxiter
    t2series = [ t2series t2 ];
	% reduce t^2(A) by, for given i,j, symmetrizing: C=AG, where G is a Givens
	% rotation such that c_ij = c_ji, then applying givens rotation such that
	% B=ZAR=ZAGZ' has b_ij=b_ji=0
	
    % find i,j
    [ i, j ] = feval( Params.fun, At );
    
	% c/s = r = (a_ii + a_jj)/(a_ji-a_ij)
	if At(i,j) ~= At(j,i) % FIX FIX should take small values into account here....
        r = ( At(i,i) + At(j,j) ) / ( At(j,i) - At(i,j) );
        s = sqrt( 1/(1+r^2) );
        c = r*s;
        % G = [ ... c .. s ...; .... -s .. c ..];
        % C = A*G
        temp   = At(:,i);
        At(:,i) =  c*temp - s*At(:,j);
        At(:,j) =  s*temp + c*At(:,j);
	else
        c = 1;
        s = 0;
	end;
	%
	if At(i,j) ~= 0  % FIX FIX should take small values into account here....
        % Z' = [ c s; -s c ]   (Z transpose!)
        % FIX FIX check: C_ii - C_jj large enough,  tayler expand.....
        tau = 2*At(i,j)/(At(i,i)-At(j,j));
        t   = (sqrt( 1 + tau^2 ) - 1)/tau;
        cr  = 1/sqrt(1+t^2);
        sr  = -t*cr;
        % multiply At from right with Z' = Givens(i,j,c & s)
        temp   = At(:,i);
        At(:,i) =  cr*temp - sr*At(:,j);
        At(:,j) =  sr*temp + cr*At(:,j);
	
        % multiply B from left with Z = Givens(i,j,c & s)'
        temp   = At(i,:);
        At(i,:) =  cr*temp - sr*At(j,:);
        At(j,:) =  sr*temp + cr*At(j,:);
        At(i,j) = 0;
        At(j,i) - 0;
	end;
    t2   = norm( At, 'fro' ) - norm( diag(At),'fro');
    iter = iter + 1;
    if mod( iter, Params.dispiter ) == 0
        disp( sprintf( 'Iter %6d, t^2: %12g', iter, t2 ) );
    end;
end;
if mod( iter, Params.dispiter ) ~= 0
    disp( sprintf( 'Iter %6d, t^2: %12g', iter, t2 ) );
end;

t2series(1) = [];   % delete initial dummy element 
svals       = sort( abs( diag( At ) ) );


% different functions to obtain i,j, passed in by name
function [ i,j ] = maxabs( At )
    Afind = At - diag(diag(At));
    [i,j] = find( abs(Afind) == max( max( abs( Afind ) ) ) );
    i=i( 1 );
    j=j( 1 );

function [ i,j ] = nextij( At )
	persistent II JJ
	if ~size(II)
        II = 1;
        JJ = 2;
	end;
   % Afind = At - diag(diag(At));
   % while At(II,JJ)^2 + At(JJ,II)^2 < 2*norm( Afind, 'fro' )/( ( size( At, 2 ) )^2 ) 
        JJ = JJ + 1;
        if JJ >= size( At, 2 )
            II = II + 1;
           if II >= size( At, 2 )
               II = 1;
            end;
            JJ = II + 1;
        end;
    % end;
	i = II;
    j = JJ;
    
    
function [ i,j ] = maxsumsq( At )
i = 2; j = 1; m = -1; n = size( At, 2 );
for ii = 2:n
    for jj = 1:(ii-1)
        if At(ii,jj)^2 + At(jj,ii)^2 > m
            i = ii;
            j = jj;
            m = At(i,j)^2 + At(j,i)^2;
        end
    end
end

function [ i,j ] = maxabsprod( At )
i = 2; j = 1; m = -1; n = size( At, 2 );
for ii = 2:n
    for jj = 1:(ii-1)
        if abs( At(ii,jj) * At(jj,ii) ) > m
            i = ii;
            j = jj;
            m = abs( At(ii,jj) * At(jj,ii) );
        end
    end
end

    
function [ i,j ] = maxsumabs( At )
i = 2; j = 1; m = -1; n = size( At, 2 );
for ii = 2:n
    for jj = 1:(ii-1)
        if abs( At(ii,jj) ) + abs( At(jj,ii) ) > m
            i = ii;
            j = jj;
            m = abs( At(i,j) ) + abs( At(j,i) );
        end
    end
end