function [ v, beta, mu ] = house( x )
% HOUSE computes v with v(1) = 1 and beta and mu such 
% that H = I - beta v v' (symmetric, orthogonal) reflects x onto  
% (the positive part of) the first coordinate axis: H*x = mu*e_1,
% with mu = 2-norm( x ) = sqrt( x'x ) >= 0, and v = x - mu * e_1. 
%
% $Id$

sigma = x(2:end)'*x(2:end);                          % flops: 2N
v     = x;                                           % mem copy: N
if isempty( sigma ) | sigma == 0
    beta = 0;
    v(1) = 1;
    mu   = abs( x(1) );   % == norm( x );
else    % note: here, we always choose v = x - norm(x) * e1
    mu = sqrt( x(1)^2 + sigma );  % == norm( x )
    if x(1) <= 0
        v(1) = x(1) - mu;
    else
        v(1) = -sigma/( x(1) + mu );  % this is also x(1) - mu,
        % but because of potential cancellation written as  
        % (x(1)^2 - mu^2)/(x(1)+mu)
    end;
    beta = 2*v(1)^2/( sigma + v(1)^2 );             % flops: 5
    v = v./ v(1);                                   % flops: N
end;

% HH reflector now: I - beta * q * q'
% flop count: about 3N (+ 10)