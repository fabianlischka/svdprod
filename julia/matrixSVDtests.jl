# matrixSVDtests.jl 
# 2013-10-08, Fabian R Lischka, at HackerSchool, New York, NY
# This file contains some test functions for the functions provided in matrixSVD.jl, 
# eg house, symtridhh, bidiaghh, etc.



# require("../src/matrixSVD.jl")
# using matrixSVD


function logErr(msg)
    println("SVDTestErr: ", msg)
end



# These genXYZ functions generates MxN testing matrices with certain structure and features,
# possibly random.
# Note: they all take M,N, with M >= N, for more consistent calling conventions, 
# but M or N might be ignored (and a square matrix returned)

function genOrthogonal(M)
    return full(qrfact!(rand(M,M))[:Q])
    # FIXFIXFIX - 
# Note: should create functions "randomOrthoPreMult", "randomOrthoPostMult"
# http://blogs.sas.com/content/iml/2012/03/28/generating-a-random-orthogonal-matrix/
# Note: see http://www.math.wsu.edu/faculty/genz/papers/rndorth.ps
end

function genUpperOnes(M,N)
    expl = "ones on the diagonal, -1 above. Eigenvals are 1, only 1 eigenvect, condition number quite bad."
    return triu(-ones(M,M)) + 2*diagm(ones(M))
end

function genNegDiag(M,N)
    expl = "upper triangular, large negative diagonal. Reasonably well conditioned."
    return triu( 10*rand( M,M ) - 1 ) - 20 * diagm( ones( M ) )
end

function genUniform(M,N)
    expl = "elements uniformly random -1..1. Complex evals."
    return 2*rand( M, N ) - 1
end

function genExtSize(M,N)
    expl = "elements of very different size. Condition ok"
    return exp(16*rand( M, N )-15).*(rand( M, N )-.5)
end

function genExtSV(M,N)
    expl = "very different singular values: 1/2 .. 1/(2^N)"
    return genOrthogonal(N) * diagm(2.0.^[-1:-1:-N]) * genOrthogonal(N)
end



function genHilb(N)
    expl = "Hilbert matrix: example of a poorly conditioned matrix. H(i,j) = 1/(i+j-1)"
    return [1.0/(i+j-1) for i=[1:N], j=[1:N]]
end

matrixGenerators = [
    (genUpperOnes, []),
    (genHilb, [])
]

# tests whether a matrix Q is orthogonal (i.e. Q'*Q = I within machine tolerance)
function testOrthogonality(Q)   
    # fnorm == sqrt(sum((Q*Q'-eye(size(Q,1))).^2))
    fnorm = normfro(Q*Q'-eye(size(Q,1)))
    cutoff = prod(size(Q))*eps()
    if fnorm > cutoff
        logErr("testSymtridhh: Frobenius norm too big: $fnorm, or relative $(fnorm/cutoff).")
        return false
    end
    return true
    # or useBLAS asum?
end


function testHouse( x ) # where x is a vector
    # call (v,beta,mu) = house(x)
    # build Householder matrix H = I - beta v v'
    # want to test: 
    # H symmetric
    # H orthogonal
    # y = H x = mu e_1 (ie y[1] == mu, y[k] == 0 for k > 1 to machine tolerance), 
    # and y[1] >= 0 and y[1] == two-norm of x
    (v,beta,mu) = house(x)
    N = size(x,1)
    # non-BLAS
    H = eye(size(x)[1]) - beta .* (v * v')  # useBLAS?

    # BLAS
    HH = eye(size(x)[1])
    Base.LinAlg.BLAS.syrk!('U','N',-beta,v,1.0,HH)
    symmetrize!(HH)

    if !issym(H)
        logErr("testHouse: H not symmetric")
        return false
    end
    if !testOrthogonality(H)
        logErr("testHouse: H not orthogonal")
        return false
    end
    y = H*x # useBLAS
    x2norm = norm(x)
    if abs(y[1] - mu) > N*eps()
        logErr("testHouse: mu not equal H*x [1]")
        return false
    end
    if abs(x2norm - mu) > N*eps()
        logErr("testHouse: mu not equal two-norm of x")
        return false
    end
    if norm(y[2:]) > N*eps()
        logErr("testHouse: H*x [2:] not sufficiently close to zero")
        return false
    end
    return true
end


function testSymtridhh( A )
    # tests the tridiagonalization of the symmetric matrix A using 
    # Q = symtridhh!(A; requireQ = true)  where A is overwritten with T
    # In particular, tests that:
    # Q*T*Q' ≈ A (the frobenius norm of the difference is below N*N*eps)
    # Q orthogonal
    # T symmetric, and tridiagonal

    T = copy(A)
    Q = symtridhh!(T; requireQ = true)
    fnorm = normfro(Q*T*Q'-A)
    cutoff = prod(size(A))*eps()
    if fnorm > cutoff
        logErr("testSymtridhh: Q*T*Q' != A, Frobenius norm too big: $fnorm, or relative $(fnorm/cutoff).")
        return false
    end
    if !testOrthogonality(Q)
        logErr("testSymtridhh: Q not orthogonal")
        return false
    end
    if !issym(T)
        logErr("testSymtridhh: T not symmetric")
        return false
    end
    # test that T is tridiagonal by setting the tri diagonals to zero and check for zero
    for k in [-1,0,1]
        T[diagind(T,k)] = 0.0
    end
    fnorm = normfro(T)
    if fnorm > cutoff
        logErr("testSymtridhh: T not tridiagonal, Frobenius norm too big: $fnorm, or relative $(fnorm/cutoff).")
        return false
    end
    return true
end



function testBidighh(A)
    # tests the bi-diagonalization of the matrix A (with M >= N) using 
    # (U,V) = bidighh!( A; requireUV = true )  where A is overwritten with B
    # In particular, tests that:
    # U*B*V' ≈ A (the frobenius norm of the difference is below N*N*eps)
    # U and V orthogonal, size U = MxM, size V = NxN
    # B upper bidiagonal

    B = copy(A)
    (U,V) = bidighh!(B; requireUV = true)
    fnorm = normfro(U*B*V'-A)
    cutoff = prod(size(A))*eps()
    if fnorm > cutoff
        logErr("testBidighh: U*B*V' != A, Frobenius norm too big: $fnorm, or relative $(fnorm/cutoff).")
        return false
    end
    if !testOrthogonality(U)
        logErr("testBidighh: U not orthogonal")
        return false
    end
    if !testOrthogonality(V)
        logErr("testBidighh: V not orthogonal")
        return false
    end
    # test that B is upper bidiagonal by setting the bi diagonals to zero and check for zero norm
    for k in [0,1]
        B[diagind(B,k)] = 0.0
    end
    fnorm = normfro(B)
    if fnorm > cutoff
        logErr("testBidighh: B not bidiagonal, Frobenius norm too big: $fnorm, or relative $(fnorm/cutoff).")
        return false
    end
    return true
end
