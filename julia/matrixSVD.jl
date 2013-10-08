# matrixSVD
# 2013-10-08, Fabian R Lischka, at HackerSchool, New York, NY
# this file contains some functions for matrix manipulation,
# mostly orthogonal transformations, and the SVD
# NOTE: this is not production quality code, use the built in functions


# function house( x :: Vector{Number}) 
# takes x : Vector, returns (v,beta,mu), v same type as x, beta, mu scalars
# HOUSE allows you to compute a Householder reflection for given 
# column vector x, such that H x lies on the first coordinate axis. 
# In particular, it computes a vector v (normalized with v(1) = 1)
# and scalars beta and mu such that
# H = I - beta v v' reflects x onto (the positive part of) 
# the first coordinate axis: H x = mu e_1, where 
# e_1 = (1 0 0 0 ...) (as a column vector) and 
# mu = two norm( x ) = sqrt( x'x ) >= 0, and 
# v = x - mu * e_1. 
function house( x ) 
    sigma = dot(x[2:],x[2:]) # useBLAS                      # flops: 2N
    # FIXFIXFIX overflow! http://blogs.mathworks.com/cleve/2012/07/30/pythagorean-addition/
    # http://degiorgi.math.hr/~singer/aaa_sem/Float_Norm/p15-blue.pdf
    v     = copy(x)                                                       # mem copy: N
    if sigma == 0  # the vector x has only length one, or all elements after the first are zero
        if x[1] >= 0
            beta = 0
        else
            # Note: we make this choice so that the invariance H x = mu e_1
            # always holds (so in this case, H x = -x = abs(x(1)) e_1)
            # otherwise we get sign errors in certain cases in routines that
            # use house and rely on that behaviour (eg. if we don't multiply
            # the first column, but just set it [mu 0 0 0 0]')
            beta = 2 
        end
        v[1] = 1 
        mu   = abs( x[1] )   # == norm( x );
    else    #  note: here, we always choose v = x - norm(x) e_1
        mu = sqrt( x[1]^2 + sigma );  # == norm( x )
        if x[1] <= 0
            v[1] = x[1] - mu
        else
            v[1] = -sigma/( x[1] + mu )
            # this is also x(1) - mu, but because of potential cancellation written as  
            # (x(1)^2 - mu^2)/(x(1)+mu)
        end
        beta = 2*v[1]^2/( sigma + v[1]^2 )                                # flops: 5
        v = v ./ v[1]                                                     # flops: N
    end

    return (v,beta,mu)
end

# HH reflector now: H = I - beta .* (v * v'), and H x = mu e_1
# flop count: about 3N (+ 10)



function symtridhh!( A; requireQ = false) 
# SYMTRIDHH brings a symmetric real matrix A in similar tridiagonal form T
# using similar householder transformations, so that T = Q'AQ with Q orthogonal.
# A must be symmetric, real. Modifies A to tridiagonal T, and returns Q if desired. 
# Note: computing Q is costly, should be done only if required, by setting requireQ = true.
# calling example: 
#   Q = symtridhh!(A; requireQ = true)
#   symtridhh!(A)
# Afterwards, Q*Q' = Q'*Q = I, and QTQ' = A.

# reference: Golub, Van Loan; 3rd ed; ch. 8.3.1 and 5.1.6.

    if !issym(A) 
        error("Input A must be symmetric!")
    end

    N = size( A, 1 )

    for k = 1:(N-2)
        (v, beta, mu ) = house( A[ (k+1):N, k ] )                   # flops: 3(N-k)
        p = A[ (k+1):N, (k+1):N ] * (beta*v)   # useBLAS                     # flops: 2(N-k)^2
        w = p - beta/2.0*dot(p,v)*v                                   # flops: 4(N-k)
        A[ k, (k+2):N ] = 0                                         # mem:   N-k
        A[ k, k+1 ] = mu                                     
        if requireQ 
            # collect the v
            A[ k+1:N, k ]   = v          
            A[ k+1, k ] = beta 
        else
            # else leave A symmetric, tridiagonal
            A[ (k+2):N, k ] = 0                                     # mem:   N-k   
            A[ k+1, k ] = mu          
        end
        # FIXFIX: take advantage of symmetry here
        A[ (k+1):N, (k+1):N ] = A[ (k+1):N, (k+1):N ] - v*w' - w*v'  # FIXFIX use syrk! ? # useBLAS
        # flops: without taking advantage of sym:                   # flops: 3(N-k)^2
    end

    if requireQ 
        # accumulate Q, backwards (smaller flop count)
        Q = eye( N )
        for k = (N-2):-1:1
            v           = A[ k+1:N, k ]
            beta        = v[ 1 ]
            v[ 1 ]      = 1
            Q[k+1:N, k+1:N] = Q[ k+1:N, k+1:N ] - beta .* (v * v') * Q[ k+1:N, k+1:N ]
            # and fix the A (eliminate the accumulated vs, be symmetric)
            A[k+1:N, k] = A[k, k+1:N ]'
        end
        return Q
    end
    # else return nothing
    # A is modified

# total flops: sum k=1:N of 5(N-k)^2 = sum k=1:N of 5N^2
# total flops: = 5/3 N^3 + lower order terms
end



function bidighh!( A; requireUV = false )
# BIDIGHH computes bidiagonal B = U'AV using householder reflections
# suppose A MxN, M >= N, A real. This routine finds orthogonal U MxM, V NxN 
# such that B = U'AV is bidiagonal (the diagonal and one superdiagonal).
# A is overwritten by B. If requireUV is false, then U, V are not explicitly formed.

# reference: Golub, Van Loan; 3rd ed; 5.4.3

    M = size( A, 1 )   # rows
    N = size( A, 2 )   # cols
    if M < N
        error( "bidighh: A must have M >= N (number of rows greater or equal the number of columns)" )
    end

    betas   = zeros( N, 1 )
    gammas  = zeros( N-2, 1 )
    B = A # note: not a deep copy, ie we modify input A. This assignment just for convenience
    for k = 1:N
        # determine HH reflector for column B(k:M,k)
        ( v, beta, mu ) = house( B[k:M,k] )            # flops: 3(M-k) 
        betas[ k ]      = beta
        #% B(k,k)          = mu;                           % = norm(x)
        B[k:M,k:N]      = B[k:M,k:N] - beta .* (v * v') * B[k:M,k:N] # FIXFIXFIX useBLAS
                                                       # flops: 4*(N-k)*(M-k)
        if requireUV
            B[(k+1):M,k]    = v[2:end]
        else
            # or, if you don't want to collect it, B((k+1):M,k) = 0;
            B[(k+1):M,k]    = 0
        end
        if k < N - 1
            ( v, beta, mu ) = house( B[k,(k+1):N]' )   # flops: 3(N-k)
            gammas[ k ]     = beta
            # B(k,k+1)        = mu;                        % = norm(x)
            # note that we need to multiply from the right here:
            B[k:M,(k+1):N]  = B[k:M,(k+1):N] - beta .* B[k:M,(k+1):N] * (v * v')
                                                       # flops: 4*(N-k)*(M-k)
            if requireUV
                B[k,(k+2):N]    = v[2:end]'
            # or, if you don't want to collect it, B(k,(k+2):N) = 0;
            else
                B[k,(k+2):N]    = 0
            end
        end
    end

    # total flops (without computing U, V)
    # sum k = 1:N of 8(N-k)(M-K) + 6(M-k)
    # = 4 MN^2 - 4/3 N^3 = 4 N^2 (M-N/3), or, 
    # for M=N, 8/3 N^3

    if requireUV
        # compute U and V
        v = zeros( M, 1 )
        U = eye( M )
        for k = N:-1:1
            # HH reflector now: I - beta * v * v'. Note: v(1) = 1
            v[k]        = 1.0
            v[k+1:M]    = B[(k+1):M,k]
            U[k:M,k:M]  = U[k:M,k:M] - betas[ k ] .* (v[k:M] * v[k:M]   ') * U[k:M,k:M]
            # flop count: 4*(M-k)^2
            B[k+1:M,k]  = 0.0
        end
        V = eye( N )
        for k = N-1:-1:2
            v[k]        = 1.0
            v[k+1:N]    = B[k-1,(k+1):N]'
            V[k:N,k:N]  = V[k:N,k:N] - gammas[ k-1 ] .* (v[k:N] * v[k:N]') * V[k:N,k:N]
            # flop count: 4*(N-k)^2
            B[k-1,k+1:N]= 0.0
        end
        return ( U, V )
    end
    # if you want to return B in two diag form B = [ diag(B) [ diag(B,1); 0 ] ];

end
# flop count for V: 4/3 N^3
# flop count for U: 4( NM^2 - MN^2 + N^3/3 )
# total flops (with computing U,V)
# 4 NM^2 + 4/3 N^3, or for M==N, 16/3 N^3 = 5.3 N^3


