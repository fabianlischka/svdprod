function [ v, beta, mu ] = house( x )
% HOUSE computes v with v(1) = 1 and beta such that F = I - beta v v' is
% orthogonal, and Fx = sqrt( x'x ) e_1
% in other words, F is the reflector that reflects x onto  
% (the positive part of) the first coordinate axis
% v = x - norm(x) * e1
%
% $Id$

sigma = x(2:end)'*x(2:end);                          % flops: 2N
v     = x;                                           % mem copy: N
v(1)  = 1;  % for the case sigma == 0
if sigma == 0
    beta = 0;
    mu = abs( x(1) );   % == norm( x );
else    % note: here, we always choose v = x - norm(x) * e1, ie v(1) always < 0, and Fx(1)>0
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