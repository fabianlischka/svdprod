function [Q,R] = qrhh( A )
% % QR-decomposition, Householder reflections
% $Id$

M = size( A, 1 );   % rows
N = size( A, 2 );   % cols
if M < N
    err( 'A must have number rows >= number columns' );
end;

Q = eye(M, N);      % Q will be M x N
R = A;              % R will be N x N .. later
betas = zeros( N, 1 );

for k = 1:N
    % determine HH reflector for column A(k:M,k)
    % [ v, beta ] = house( A(k:M,k) );   
    % we will do it explicitly here - flops: 3(M-k)
    v  = R(k:M,k);
    x1 = v(1);
    sigma = v(2:end)'*v(2:end);
    if isempty( sigma ) || sigma == 0
        R(k,k) = x1;  % = norm(x). No need to update rest of matrix (or betas), since beta = 0
        % NOTE: used to be abs(x1) = norm(x), but now reverted to x1. Why?
        % because if we flip sign here to maintain positive diagonal on R,
        % we would need to flip sign in Q, but there is no way to store
        % that informaion! (since beta = 0...)
        % collect Qi: v is zero
        R((k+1):M,k) = 0;
	else    % note: here, we always choose v = x - norm(x) * e1, ie v(1) always < 0
        mu = sqrt( x1^2 + sigma );  % = norm( x )
        if x1 <= 0
            v(1) = x1 - mu;
        else
            v(1) = -sigma/( x1 + mu );  % this is also x(1) - mu,
            % but because of potential cancellation written as  
            % (x(1)^2 - mu^2)/(x(1)+mu)
        end;
        betas( k ) = 2*v(1)^2/( sigma + v(1)^2 );
        v = v./ v(1);
        % HH reflector now: I - beta * v * v'. Note: v(1) = 1
        
        R(k,k) = mu;                            % = norm(x)
                                                % flops until here: 3(M-k)
        R(k:M,(k+1):N) = R(k:M,(k+1):N) - betas( k ) * v * v' * R(k:M,(k+1):N);
        % (note: if we computed (v*v')*R, we would neet 2MN^2 flops, but
        % with v*(v'*R) we only need 4NM!!)
                                                % flops: 4*(N-k)(M-k) 
        % collect Qi
        R((k+1):M,k) = v(2:end);
	end;
end;

% compute Q, and zero out R
v = zeros( M, 1 );  % note: without this, v is switched to a row vector!!
for k = (N-1):-1:1
    % HH reflector now: I - beta * v * v'. Note: v(1) = 1
    v(k)=1;
    v((k+1):M) = R((k+1):M,k);
    Q(k:M,k:N) = Q(k:M,k:N) - betas( k ) * v(k:M) * transpose(v(k:M)) * Q(k:M,k:N);
    R((k+1):M,k) = 0;
end;
R = R(1:N,:);

% flop count: without computing Q: 
% sum k=1:N of 4(N-k)(M-k) = 4 sum k=1:N k(M-N+k)
% = 4 sum k^2 + 4(M-N) sum k = 4/3 N^3 + 4/2(M-N)N^2
% = 2 N^2(M-N/3)
