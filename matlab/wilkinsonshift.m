function mu = wilkinsonshift( a1, b, a2 )
% WILKINSONSHIFT returns the eigenvalue of the symmetric matrix 
% (a1, b; b, a2) that is closer to a2. Some care is taken to avoid
% cancellation.
% TS N-1,1    a1
% TS N  ,1    a2
% TS N-1,2    b

% reference: Golub, Van Loan; 3rd ed; ch. 8.3.4
% $Id$

d       = ( a1 - a2 ) / 2;    % flops: 2
% mu      = a2 + d - sign( d ) * sqrt( d^2 + b^2 );
% note that ( d - sign(d)*sqrt(d^2+b^2) ) * ( d + sign(d)*sqrt(d^2+b^2) ) =
%           = d^2 - (d^2 + b^2) = - b^2, so we can avoid cancelation
%           between d and the square root, particularly when b^2 (which
%           is the last superdiagonal element) is small - which we hope!
if d == 0
    % then a1==a2 =: a, evals are a +- abs(b) 
    % we choose absolutely larger to avoid cancellation
    if a2 > 0
        mu      = a2 + abs( b );
    else
        mu      = a2 - abs( b );
    end
else
    mu      = a2 -  b^2 / ( d + sign(d)*sqrt( d^2 + b^2 ) );
end;