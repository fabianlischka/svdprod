% % QR-decomposition, Householder reflections
% $Header$

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
K = min( M, N );    % min = steps

Q = zeros(M, K );   % Q will be M x K
R = A(1:K,:);       % R will be K x N

for k = 1:(K-1)
    % determine HH reflector for column A(k:K,k)
    % [ v, beta ] = house( A(k:K,k) );
    v  = R(k:K,k);
    x1 = v(1);
    sigma = v(2:end)'*v(2:end);
    if sigma == 0
        R(k,k) = abs( x1 );  % = norm(x). No need to update rest of matrix, since beta = 0
	else    % note: here, we always choose v = x - norm(x) * e1, ie v(1) always < 0
        mu = sqrt( x1^2 + sigma );  % = norm( x )
        if x1 <= 0
            v(1) = x1 - mu;
        else
            v(1) = -sigma/( x1 + mu );  % this is also x(1) - mu,
            % but because of potential cancellation written as  
            % (x(1)^2 - mu^2)/(x(1)+mu)
        end;
        beta = 2*v(1)^2/( sigma + v(1)^2 );
        v = v./ v(1);
        % HH reflector now: I - beta * v * v'. Note: v(1) = 1
        
        R(k,k) = mu;                            % = norm(x)
        R(k:K,(k+1):K) = R(k:K,(k+1):K) - beta * v * v' * R(k:K,(k+1):K);
        R((k+1):K,k) = 0;
	end;
end;