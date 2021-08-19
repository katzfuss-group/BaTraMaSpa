using Random
using StatsBase
using Plots
using JLD
using LinearAlgebra
using KernelDensity
include("./densities.jl")
include("./chains.jl")
include("./reference_measures.jl")
include("./samplers.jl")
# defining the target distribution
q = 6
n = 2 ^ q
μ = 0.05
ν = 1.0
T = 25

# defining the sampler
δ = 0.25
k_max = 10000000

# decreasing weights according to temperature
μ /= T
ν /= T
lin_inds = LinearIndices((n, n))

# creates sampler for density 
sampler = PCNMHSampler(create_FD_laplacian_reference_measure(q, μ),
                       CahnHilliardDensity(ν),
                       δ)


@time chain, acceptance_probability = sample_chain(sampler, k_max, round(Int, k_max / 10))

@show acceptance_probability

save("out/unconditional_chain_kmax_$(k_max)_nu_$(ν)_mu_$(μ)_T_$(T)_q_$(q)_delta_$(δ).jld", 
     "chain", chain,
     "μ", μ,
     "ν", ν,
     "T", T,
     "k_max", k_max,
     "q", q,
     "δ", δ)