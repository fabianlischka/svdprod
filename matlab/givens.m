function [ c, s ] = givens( a, b )
% GIVENS computes c and s such that [ c -s; s c ] [ a b ]' = [ r 0 ]'
% protected against overflow

% reference: Golub, Van Loan; 3rd ed; 5.1.8
% $Id$

if b == 0
    c = 1; s = 0;
else
    if abs( b ) > abs( a )
        tau = -a/b; s = 1/sqrt( 1+tau^2 ); c = s*tau;
    else
        tau = -b/a; c = 1/sqrt( 1+tau^2 ); s = c*tau;
    end;
end;

% flop count: 5 and 1 sqrt