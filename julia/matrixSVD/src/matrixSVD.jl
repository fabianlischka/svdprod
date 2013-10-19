# matrixSVD
# 2013-10-08, Fabian R Lischka, at HackerSchool, New York, NY
# this file contains some functions for matrix manipulation,
# mostly orthogonal transformations, and the SVD
# NOTE: this is not production quality code, use the built in functions

# main SVD routines
export gksvd, gksvdprod
# helpers (might be useful)
export bidighh!, bidighh, house, wilkinsonshift, givens, bidiagshortform
# others
export symtridhh!, symtridhh


# internal funcs: gksvdstep!, gksvdsteps!


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
function house( x :: Vector) # FIXFIXFIX - type?
    sqrtSigma = norm(x[2:])
    v     = copy(x)                                 # mem copy: N
    if sqrtSigma == 0  # the vector x has only length one, or all elements after the first are zero
        if x[1] >= 0
            β = 0.0
        else
            # Note: we make this choice so that the invariance H x = μ e_1
            # always holds (so in this case, H x = -x = abs(x(1)) e_1)
            # otherwise we get sign errors in certain cases in routines that
            # use house and rely on that behaviour (eg. if we don't multiply
            # the first column, but just set it [μ 0 0 0 0]')
            β = 2.0
        end
        v[1] = 1.0
        μ   = abs( x[1] )   # == norm( x );
    else    #  note: here, we always choose v = x - norm(x) e_1
        μ = hypot(x[1],sqrtSigma)  # was: sqrt( x[1]^2 + sigma )  # == norm( x )
        if x[1] <= 0
            v[1] = x[1] - μ
        else
            v[1] = -sqrtSigma^2/( x[1] + μ )    # was:  -sigma/( x[1] + μ )
            # this is also x(1) - μ, but because of potential cancellation written as
            # (x(1)^2 - μ^2)/(x(1)+μ)
        end
        β = 2*v[1]^2/( sqrtSigma^2 + v[1]^2 )       # flops: 5ish FIXFIXFIX: overflow
        v = v ./ v[1]                               # flops: N
    end
    return (v,β,μ)
end

# HH reflector now: H = I - β .* (v * v'), and H x = μ e_1
# flop count: about 3N (+ 10)




function bidighh!( A :: Matrix; requireUV = false )
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
            U[k:M,k:M]  = U[k:M,k:M] - betas[ k ] .* (v[k:M] * v[k:M]') * U[k:M,k:M]
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


bidighh( A :: Matrix; requireUV = false ) = bidighh!( copy(A); requireUV = requireUV )


function bidigprod( A; requireUV = false )
# BIDIGHH computes bidiagonal B=U' AAA V using householder reflections.
# suppose A NxNxK, real. Without ever computing the product
# AAA = A(:,:,K)*...*A(:,:,2)*A(:,:,1), this routine finds
# orthogonal U NxN, V NxN such that B = U' AAA V is bidiagonal
# (the diagonal and one superdiagonal).
# B is returned in short format, ie the diagonal
# is stored in the first col of B, the superdiagonal in the
# second col of B, and B(end,2) == 0

# reference: Golub, Solna, van Dooren: Computing the SVD
# of a General Matrix Product/Quotient,
# SIAM J. MATRIX ANAL. APPL., Vol. 22, No. 1, pp 1-19
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.51.4387&rep=rep1&type=pdf

    N   = size( A, 1 )
    if size( A, 2 ) ~= N
        error( "bidigprod: matrices in A must be square (ie size(A,1)=size(A,2))" )
    end

    K   = size( A, 3 )
    if requireUV
        V   = eye( N )
        A[:,:,K+1] = V    # this is U - FIXFIXFIX why?
    end;

    for t=1:N-1
        for k=1:K
            ( v, beta, mu ) = house( A[ t:N, t, k ] )   # flops: 3(N-t)
            # I-beta vv' is Q^t_k, applied to t:N
            # multiply Ak from left...
            A[ t:N,t:N, k ] = A[ t:N,t:N, k ] - beta.* v*v' * A[ t:N,t:N, k ]
                                                        # flops: 4K(N-t)^2
            # FIXFIXFIX: look into these householder updates, want to maintain symmetry, but be fast
            if k<K || requireUV    # execute almost always
                # ...and Ak+1 from right - here we need to multiply all columns,
                # since we don't have the zeros needed in here
                A[ :,t:N, k+1 ] = A[ :,t:N, k+1 ] - A[ :,t:N, k+1 ]*v*v'.*beta;
                                                        # flops: 4(K-1)N(t-N)
                                                        # with vector: +4N(t-N)
            end
        end

        if t < N-1
            # now, determine the t-th row of the product
            row = A[ t, t:N, K ]
            for k = K-1:-1:1
                row = row * A[ t:N, t:N, k ]            # (flops: 2(N-t)^2)
            end
                                                        # flops: 2K(N-t)^2
            # and determine householder to eliminate it!
            # (the diagonal element remains untouched (and is thrown away))
            ( v, beta, mu ) = house( row[ 1, 2:end ]' ) # flops: 3(N-t)
            # multiply A( 1 ) and V from the right
            # need to multiply full columns (all rows), since might be filled
            A[ :, t+1:N, 1 ] = A[ :,t+1:N, 1] - A[ :,t+1:N, 1 ]*v*v'*beta
                                                        # flops: 4N(t-N)
            if requireUV
                V[ :,t+1:N ] = V[ :,t+1:N ]   - V[ :, t+1:N ]  *v*v'*beta
                                                        # flops: 4N(t-N)
            end
        end
        # total flops per t:  with vectors (8+4K) N(t-N) + 6K (N-t)^2
        #                     = (4+4K) N^3
        #                     without vectors: (4K) N^3
    end
    # now B = A(:,:,K+1)'*AA*V = A(:,:,K)*...*A(:,:,2)*A(:,:,1)

    q   = ones(N,1)      # this will be the diagonal of the product B
    e   = zeros(N-1,1)   # superdiagonal of the product B
    # note: to compute q, e, we only need diagonal and superdiag of A !
    for k=1:K
        d   = diag( A[:,:,k] )
        e   = e .* d(1:end-1) + q(2:N) .* diag( A[:,:,k], 1 )
        q   = q .* d
    end

    if requireUV
        U = A[:,:,K+1]
    end
    B = [ q [e; 0] ]
    return (B,U,V)
    # note: B returned in short format
end



function wilkinsonshift( a1, b, a2 )
# WILKINSONSHIFT returns the eigenvalue of the symmetric 2x2 matrix
# (a1, b; b, a2) that is closer to a2. Some care is taken to avoid
# cancellation.
# % TS N-1,1    a1
# % TS N  ,1    a2
# % TS N-1,2    b

# reference: Golub, Van Loan; 3rd ed; ch. 8.3.4

    d  = ( a1 - a2 ) / 2           # flops: 2
# % note that ( d - sign(d)*sqrt(d^2+b^2) ) * ( d + sign(d)*sqrt(d^2+b^2) ) =
# %           = d^2 - (d^2 + b^2) = - b^2, so we can avoid cancelation
# %           between d and the square root, particularly when b^2 (which
# %           is the last superdiagonal element) is small - which we hope!
    if d == 0
        # then a1==a2 =: a, evals are a +- abs(b)
        # we choose absolutely larger to avoid cancellation
        if a2 > 0
            μ = a2 + abs( b )
        else
            μ = a2 - abs( b )
        end
    else
        μ = a2 -  b^2 / ( d + sign(d)*hypot( d, b) )
    end
    return μ
end
# flop count: around 14


function givens( a, b )
# GIVENS computes c and s such that [ c s; -s c ]' [ a b ]' = [ r 0 ]'
# protected against overflow

# reference: Golub, Van Loan; 3rd ed; 5.1.8

    if b == 0
        return (1.0,0.0)
    else
        if abs( b ) > abs( a )
            tau = -a/b
            s = 1/sqrt( 1+tau^2 )
            c = s*tau
        else
            tau = -b/a
            c = 1/sqrt( 1+tau^2 )
            s = c*tau
        end
    end
    return (c,s)
end
# flop count: 5 and 1 sqrt


function bidiagshortform(B)
    # turns the upper NxN bidiagonal matrix B into "short form" Nx2
    # with the diagonal in the 1st, and the superdiagonal + [0] in the 2nd col
    return [ diag(B) [ diag(B,1); 0 ] ]
end

function gksvdstep!( B, U, V )
# Given an (upper) bidiagonal real B, NxN, with no zeros on the diagonal or
# superdiagonal, GKSVDSTEP (Golub-Kahan Singular Value Decomposition Step)
# computes B+  = U' B V, where U and V NxN are orthogonal matrices,
# and VV is the orthogonal matrix that would have been obtained by applying
# an implicit QR step with Wilkinson shift (QRIMSTEP) to T = B'B,
# ie T+ = B+'B+ = V'B'U U'BV = V' B'B V = V'TV
# and we want T+ triangular, and V*e1 proportional to T-mu*I

# B is given in short format, ie. the diagonal is in the first col B(1:N,1),
# the superdiagonal in the second col B(1:N-1, 2).

# Note: internal function, no checks on the inputs performed

# reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.1

    N = size( B, 1 )

    # determine the shift mu.
    # mu should be the eigenvalue (closer to the lower right) of
    # ( B(N-1,1)^2 + B(N-2,2)^2    ,  B(N-1,1) * B(N-1,2)    )
    # ( B(N-1,1) * B(N-1,2)        ,  B(N,1)^2 + B(N-1,2)^2  )
    if N > 2
        mu = wilkinsonshift( B[N-1,1]^2+B[N-2,2]^2,  B[N-1,1]*B[N-1,2],  B[N,1]^2+B[N-1,2]^2 )
    else if N == 2
            mu = wilkinsonshift( B[N-1,1]^2,  B[N-1,1]*B[N-1,2],  B[N,1]^2+B[N-1,2]^2 )
        else
            error("Cannot perform GK SVD step on B with N < 2")
        end
    end

    x  = B[1, 1]^2 - mu
    y  = B[1, 1] * B[1, 2]
    
    # B+ = U' B V =   U(N-1)' ... U2' U1' B G1 V2 V3 ... V(N-1), with G1*e1 prop T-mu*I
    for k = 1:N-1
        # compute Givens rotation, around k,k+1 FROM RIGHT, G = [c s; -s c]
        c, s     = givens( x, y )                           # flops: 5 + 1 sqrt
        # rotate B: want B_(1/2) = B G, V+ = V G
        G           = [ c s; -s c ]
        # update V - we update the full column, since we don't know what was
        # passed in. FIXFIXFIX
        V[:, k:k+1] = V[:, k:k+1] * G                       # flops: 6 * N
        # for B in ordinary format, we want:
        # B( k-1:K, k:k+1 ) = B( k-1:K, k:k+1 ) * G;

        if k > 1
            # for debugging we can compute the zero'ed out element:
            # zerod   = s * B(k-1,2) + c * bulge
            B[ k-1, 2 ] =  c * B[ k-1, 2 ] - s * bulge      # flops: 3
        end
        # chase the bulge...
        Bk          = B[ k, 1 ]
        bulge       = -s * B[ k+1, 1 ]                      # flops: 1
        B[ k, 1 ]   =  c * Bk - s * B[ k, 2 ]               # flops: 3
        B[ k, 2 ]   =  s * Bk + c * B[ k, 2 ]               # flops: 3
        B[ k+1, 1 ] =  c * B[ k+1, 1 ]                      # flops: 1

        # now, bulge is at k+1, k
        x           = B[ k, 1 ]
        y           = bulge       # which we want to eliminate...

        # compute Givens rotation, around k,k+1 FROM LEFT, G = [c s; -s c]
        c, s        = givens( x, y )                        # flops: 5 + 1 sqrt
        # rotate B: want B_+ = G' B_(1/2), U+ = U G
        G           = [ c s; -s c ]
        # update U (from right, since transposed)
        # we update the full column, since don't know what was passed in
        U[:, k:k+1] = U[:, k:k+1 ] * G                      # flops: 6 * K
        # for B in ordinary format, we want:
        # B = G' * B;

        # zerod       =  s * B[ k, 1 ) + c * bulge
        B[ k, 1 ]   =  c * B[ k, 1 ] - s * bulge
        Bk2         =  B[ k, 2 ]
        B[ k, 2 ]   =  c * Bk2       - s * B[ k+1, 1 ]
        bulge       = -s * B[ k+1, 2 ]
        B[ k+1, 1 ] =  s * Bk2       + c * B[ k+1, 1 ]
        B[ k+1, 2 ] =  c * B[ k+1, 2 ]
        # now, bulge is at k, k+2

        x           = B[ k, 2 ]
        y           = bulge       # which we want to eliminate...
    end
    B[ N, 2 ] = 0
end

function gksvdsteps!( B, U, V, tol )
# GKSVDSTEPS computes the SVD of a bidiagonal B, and
# overwrites U, V such that if previously B=U'AV, now D=U'AV
# this is an internal function, so no errorchecking on input
# both B and D in short format (as column vectors)

# reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.2
# $Id$

    M = size( U, 1 )   # rows
    N = size( V, 1 )   # cols

    # partition B into an unreduced part in the middle, and a diagonal part at
    # the bottom
    q = 0
    while q < N
        q   = 0
        p = N-1
        for k = N-1:-1:1
            # note: below, we need <=, not <, otherwise we get caught in an
            # infinite loop, if these elements are exactly zero!
            if abs( B[k,2] ) <= tol*( abs( B[k,1] ) + abs( B[k+1,1] ) )
                B[k,2] = 0
                if q == N-k-1   # state one
                    q = N-k
                    p = k-1
                # else          # state three
                end
            else
                if p == k       # state two
                    p = k-1
                # else          # state three
                end
            end
        end
        if q == N-1
            q=N    # and abort...
            # NOW, with D11 = BL(1:p,1:p), D22 = BL(p+1:N-q,p+1:N-q),
            # D33 = BL(N-q+1:N, N-q+1:N), D33 is diagonal, D22 unreduced
        else        # do some work
            # first, determine if any element on the diagonal is zero
            # if so, rotate it aside
            k = p+1
            smalldiag = tol * norm( B, Inf ) # FIXFIXFIX infinity norm?
            while abs( B[ k, 1 ] ) > smalldiag  && k < N-q
                k = k+1
            end
            # now, k == N-q, or abs( B[k,1) ) <= smalldiag, or both
            if abs( B[ k, 1 ] ) <= smalldiag
                B[ k, 1 ] = 0
                # now, B[ k, 1 ) approx 0, p+1 <= k <= N-q, B[N-q,2) == 0 (note:
                # even if q == 0, by the format we have chosen)
                if k < N-q
                    # to do: if k<N-q, rotate that row away with Givens rotations from
                    # left, G(k,j), j=k+1:N-q
                    bulge  = B[k,2]
                    B[k,2] = 0
                    for j=k+1:N-q
                        c,s = givens( B[j,1], bulge )
                        # need to transpose (since mult from left), and
                        # transpose again, since need to eliminate upper
                        # element
                        # rotate B from left: B+ = G B,  G=G(k,j)
                        # zerod =  c * bulge + s * B[j,1)
                        B[j,1]  = -s * bulge + c * B[j,1]
                        bulge   =  s * B[j,2]
                        B[j,2]  =  c * B[j,2]
                        # rotate U from right (since transposed): U+ = U G'
                        G       = [ c s; -s c ]
                        U[:,[k j]] = U[:,[k j]] * G'       # flops: 6*M
                    end # for j
                else # k==N-q
                    # if k==N-q, then apply Givens from the
                    # right, G(k,j), j=N-q-1:-1:p+1
                    bulge = B[N-q-1,2]
                    B[N-q-1,2] = 0
                    for j=N-q-1:-1:p+1
                        # rotate B from right, B+ = B G;  G=G(k,j)
                        c,s  = givens( B[j,1], bulge )
                        # zerod   = s * B[j,1) + c * bulge
                        B[j,1]  = c * B[j,1] - s * bulge
                        if j>p+1
                            bulge   = s * B[j-1,2]
                            B[j-1,2]= c * B[j-1,2]
                        end;
                        # rotate V from right
                        G       = [ c s; -s c ]
                        V[:,[j k]] = V[:,[j k]] * G       # flops: 6*M
                    end
                end # if k
            else    # no element on the diagonal is zero
                # do some real work - GK SVD step for p+1:N-q
                # NOTE: the sub(A,xcoords,ycoords) create a view that is not a copy but a reference
                gksvdstep!( sub(B,p+1:N-q,:), sub(U,:,p+1:N-q), sub(V,:,p+1:N-q) )
            end
        end # if q
    end # while q<N

    # now fix sign - want cols of V such that elements of D are non-negative
    D = B[:,1]
    for k = 1:N
        if D[k] < 0
            D[k]    = -D[k]
            V[:,k]  = -V[:,k]
        end
    end
    return D    # while U,V are overwritten
end


function gksvd( A, tol = 1e-14 )
# gksvd computes the SVD of A.
# NOTE: this is an example implementation only, use the Julia library SVD for production purposes!
# Given A MxN, M>=N, and tol, this computes U MxM, V NxN orthogonal and DD diagonal
# such that A=U(DD+E)V', where E is (small) error term
# (in particular, the two norm of E should be around eps x two-norm A)
# D is returned as a Nx1 column vector (of singular values),
# DD would be = [ diag(D); zeros( M-N, N ) ]   # FIXFIXFIX

# % reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.2

    M, N = size(A)
    if M < N
        error( "A must have M >= N" )
    end
    # FIXFIXFIX: check number of dimension is 2

    # bidighh! computes bidiagonal B=U'AV using householder reflections, and writes it into A
    U,V = bidighh!( A; requireUV = true )
    # gksvdsteps computes diagonal D and updates U, V such that D=U'AV
    B = bidiagshortform( A )
    D = gksvdsteps!( B, U, V, tol )
    return (U,D,V)
end

function gksvdprod( A, tol = 1e-14 )
# gksvd computes the SVD of a product of square matrices
# NOTE: this is an example implementation only, use the Julia library SVD for production purposes!
# Given A NxNxK, for AAA := A(:,:,K)*...*A(:,:,1) this computes
# U, V NxN orthogonal, and D Nx1 such that AAA=U diag(D) V'
# Note: AAA NxN is never explicitly computed!
# Note: the only difference between gksvdprod and gksvd is in the bidiagonalization

# reference: Golub, Van Loan; 3rd ed; ch. 8.6.2, alg. 8.6.2
# reference: Golub, Solna, van Dooren: Computing the SVD
# of a General Matrix Product/Quotient,
# SIAM J. MATRIX ANAL. APPL., Vol. 22, No. 1, pp 1-19
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.51.4387&rep=rep1&type=pdf

    M, N, K = size(A)
    # FIXFIXFIX: check number of dimension is 3
    if M != N
        error( "A must have M == N" )
    end

    # bidigprod! computes bidiagonal B = U' AAA V using householder reflections
    B,U,V = bidigprod( A; requireUV = true )
    # gksvdsteps computes diagonal D and updates U, V such that D=U'AV
    D = gksvdsteps!( B, U, V, tol )
    return (U,D, V)
end



################ other stuff



# function qrsymtrid( T )
# # QR decomposition for symmetric tridiagonal matrices
# # QRSYMTRID returns Q and R such that T=QR, Q is orthogonal (and upper
# # Hessenberg), and R is upper triangular (of bandwidth 2, ie diagonal and
# # two superdiagonals).

# # reference: Golub, Van Loan; 3rd ed; ch. 8.3.3

#     N = size( T, 1 )
#     if !issym(T)
#         error( "qrsymtrid: input T must be symmetric" )
#     end;

#     R   = T
#     Q   = eye( N,N )

#     # remove 0 on subdiagonal by premultiplication with orthogonal Givens
#     # matrices Q1', Q2', ... -> obtain upper triangular R = QN-1'...Q2' Q1' T
#     # so T = QR with Q = Q1 Q2 Q3 ... QN-1
#     for k=1:(N-1)
#         # now want to find Givens rotation in k, k+1 such that
#         # R( k+1, k ) = 0. This needs to be applied to row k, k+1
#         # GIVENS computes c and s such that [ c s; -s c ]' [ a b ]' = [ r 0 ]'
#         c, s            = givens( R[k,k], R[ k+1, k ] )   # flop count: 6
#         G               = [ c s; -s c ]
#         K               = min( N, k+2 )
#         R[ k:k+1, k:K ] = G' * R[ k:k+1, k:K ]            # flop count: 3*6
#         R[ k+1, k ]     = 0     # should be close to zero anyway...
#         # accumulate Q
#         Q[1:k+1, k:k+1] = Q[1:k+1, k:k+1] * G             # flop count: (k+1)*6
#     end
#     return (Q,R)
# end
# # % flop count: without Q accumulation, around 25 N
# # % with Q accumulation: sum k=1:N of 30 + 6 k
# # % total flop count 3 N^2 + 33 N


function symtridhh!( A :: Matrix; requireQ = false)
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

symtridhh( A; requireQ = false ) = symtridhh!(copy(A); requireQ = requireQ)
