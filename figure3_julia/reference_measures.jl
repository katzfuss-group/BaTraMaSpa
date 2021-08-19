using Base: Float64
using FFTW: fft, fft!, ifft, ifft!
using SparseArrays: SparseVector, sparsevec, spdiagm, sparse
using LinearAlgebra: cholesky, ldiv!, mul!
abstract type AbstractReferenceMeasure end
# an abstract type for processes that are sampling using the FFT
abstract type AbstractFFTReferenceMeasure<:AbstractReferenceMeasure end

struct FFTReferenceMeasure<:AbstractFFTReferenceMeasure
    fft_factor::Matrix{Float64} 
    temp::Matrix{ComplexF64}
    out::Matrix{Float64}
end

function create_FFT_laplacian_reference_measure(n, μ)
    # create a function that samples n Gaussians
    h = 1 / n
    stencil = zeros(n, n)      
    # creating centered (negative) Laplacian stencil
    stencil[1, 1] = 4; stencil[1, 2] = -1; stencil[2, 1] = -1; stencil[end, 1] = -1; stencil[1, end] = -1
    # divide by squared mesh constant 
    stencil ./= h^2
    # multiply with weight of regularizer
    stencil .*= μ
    fft_factor = 1 ./ sqrt.(fft(stencil))  
    fft_factor[1, 1] = 0.0
    return FFTReferenceMeasure(real.(fft_factor), 
                               zeros(ComplexF64, size(fft_factor)),
                               zeros(size(fft_factor)))
end

function sample(measure::FFTReferenceMeasure)
    prop = fft(randn(size(measure.fft_factor)))
    prop .*= measure.fft_factor
    return vec(real.(ifft(prop)))
end 
    
function sample!(measure::FFTReferenceMeasure)
    out = measure.out; temp = measure.temp; fft_factor=measure.fft_factor
    randn!(out)
    temp .= out
    fft!(temp)
    temp .*= fft_factor
    ifft!(temp)
    @. out = real(temp)
    return vec(out)
end

# in-place application of a matrix-vector multiply with square root of covariance matrix 
function multiply_with_sqrt_C!(measure::AbstractFFTReferenceMeasure, in)
    out = measure.out; temp = measure.temp; fft_factor=measure.fft_factor
    vec(temp) .= in
    fft!(temp)
    temp .*= fft_factor
    ifft!(temp)
    @. out = real(temp)
    return vec(out)
end

function multiply_with_C!(measure::AbstractFFTReferenceMeasure, in)
    out = measure.out; temp = measure.temp; fft_factor=measure.fft_factor
    vec(temp) .= in
    fft!(temp)
    temp .*= fft_factor .^ 2
    ifft!(temp)
    @. out = real(temp)
    return vec(out)
end



# creates the sparse orthonormal matrix representing the haar wavelets
function create_haar_wavelets(q)
    n = 2 ^ q
    N = n ^ 2

    # each sparse vector is a scaling function (S) or orthogonal wavelet (W)
    S = Vector{SparseVector{Float64,Int}}[]
    W = Vector{SparseVector{Float64,Int}}[]
    # add the finest partition to S
    push!(S, SparseVector{Float64,Int}[])
    for i = 1 : N
        push!(S[end], sparsevec([i], [1.0], N))
    end 

    for k = (q - 1) : -1 : 0
        nk = 2 ^ k
        Nk = nk^2
        lin_inds = LinearIndices((nk * 2, nk * 2))
        # adding a new level to S and W
        push!(S, SparseVector{Float64,Int}[])
        push!(W, SparseVector{Float64,Int}[])
        # go over all the support regions on the next coarser scale
        for (i, j) in Tuple.(CartesianIndices((nk, nk)))
            # obtaining the indices on the next finer level that correspond to basis
            # functions on the next coarser level
            ind11 = lin_inds[(i - 1) * 2 + 1, (j - 1) * 2 + 1]
            ind21 = lin_inds[(i - 1) * 2 + 2, (j - 1) * 2 + 1]
            ind12 = lin_inds[(i - 1) * 2 + 1, (j - 1) * 2 + 2]
            ind22 = lin_inds[(i - 1) * 2 + 2, (j - 1) * 2 + 2]
            # Adding the single basis function to S
            push!(S[end], S[end-1][ind11] + 
                          S[end-1][ind21] + 
                          S[end-1][ind12] +
                          S[end-1][ind22]) 

            # Adding the three basis functions to W
            push!(W[end], S[end-1][ind11] + 
                          S[end-1][ind21] - 
                          S[end-1][ind12] -
                          S[end-1][ind22]) 

            push!(W[end], S[end-1][ind11] - 
                          S[end-1][ind21] + 
                          S[end-1][ind12] -
                          S[end-1][ind22]) 

            push!(W[end], S[end-1][ind11] - 
                          S[end-1][ind21] - 
                          S[end-1][ind12] +
                          S[end-1][ind22]) 
        end
    end 
    push!(W, SparseVector{Float64,Int}[])
    push!(W[end], only(S[end]))

    for level in W
        for basis_function in level
            basis_function .= basis_function / norm(basis_function)
        end
    end
    return reduce(hcat, reduce(vcat, W))
end

struct CholeskyReferenceMeasure
    # An orthogonal matrix providing the ordering / basis transform
    H
    # The cholesky factorization of the matrix corresponding to the stochastic dofs
    chol_A22
    # The matrix coupling deterministic and stochastic dofs 
    A21
    # The values for the last length(prescribed_values) coefficients
    prescribed_values
    # scratch space
    transformed_prop::Vector{Float64}
    prop::Vector{Float64}
    temp::Vector{Float64}
end

# μ is 𝐧𝐨𝐭 the mean, but the weight with which the laplacian energy is multiplied.
function create_FD_laplacian_reference_measure(q, μ, prescribed_values=[0.0])

    n = 2 ^ q
    N = n ^ 2  
    h = 1 / n
    deterministic_indices = (N - length(prescribed_values) + 1) : N
    stochastic_indices = 1 : (N - length(prescribed_values))

    ##################################################################
    # Construct the Laplace operator
    ##################################################################
    # Create the sparsity pattern of the Laplace operator 
    lin_inds = LinearIndices((n, n))
    row_inds = Int[]
    col_inds = Int[]
    S = Float64[]
    for i in 1 : n, j in 1 : n
        # adding self-interaction 2
        push!(col_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1, n) + 1]) 
        push!(row_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1, n) + 1]) 
        push!(S, 4.0) 
        # interaction to next element in i direction
        push!(col_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1, n) + 1]) 
        push!(row_inds, lin_inds[mod(i - 1 + 1, n) + 1, mod(j - 1, n) + 1]) 
        push!(S, -1.0)

        # interaction to previous element in i direction
        push!(col_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1, n) + 1])
        push!(row_inds, lin_inds[mod(i - 1 - 1, n) + 1, mod(j - 1, n) + 1]) 
        push!(S, -1.0)

        # interaction to next element in j direction
        push!(col_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1, n) + 1]) 
        push!(row_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1 + 1, n) + 1]) 
        push!(S, -1.0)

        # interaction to previous element in j direction
        push!(col_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1, n) + 1]) 
        push!(row_inds, lin_inds[mod(i - 1, n) + 1, mod(j - 1 - 1, n) + 1]) 
        push!(S, -1.0)
    end
    A = sparse(row_inds, col_inds, S) * μ / (h^2)

    # creates the matrix that has the (L²-normalized) Haar wavelets as columns
    H = create_haar_wavelets(q)

    chol_A22 = cholesky(Hermitian(H[:, stochastic_indices]' * A * H[:, stochastic_indices]))
    A21 = H[:, stochastic_indices]' * A * H[:, deterministic_indices]

    return CholeskyReferenceMeasure(H, chol_A22, A21, prescribed_values, zeros(N), zeros(N), zeros(length(stochastic_indices)))
end
    
function sample!(measure::CholeskyReferenceMeasure)
    prop = measure.prop; temp = measure.temp; transformed_prop = measure.transformed_prop
    N = size(measure.A21, 1) + size(measure.A21, 2)
    deterministic_indices = (N - size(measure.A21, 2) + 1) : N
    stochastic_indices = 1 : (N - size(measure.A21, 2))
    # fill the stochastic part with Gaussian randomness
    randn!(temp)
    # Multiply with inverse square root to create the (centered) conditional process
    # ldiv!(measure.chol_A22.L, temp)

    temp .= measure.chol_A22.PtL' \ temp

    @views prop[stochastic_indices] .= temp
    # writing A21 v into temp, for v the prescribed values
    mul!(temp, measure.A21, measure.prescribed_values)
    # computing the expectation
    temp .= - (measure.chol_A22 \ temp)
    # adding it to the proposal
    view(prop, stochastic_indices) .+= temp
    view(prop, deterministic_indices) .= measure.prescribed_values
    mul!(transformed_prop, measure.H, prop)
    return transformed_prop
end