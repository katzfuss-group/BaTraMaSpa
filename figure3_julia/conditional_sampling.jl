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

# the depth of the conditioning
p = 5

# the parameters for the unconditional problem
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

# loading the initial data from the unconditional problem
ld = load("out/unconditional_chain_kmax_$(k_max)_nu_$(ν)_mu_$(μ)_T_$(T)_q_$(q)_delta_$(δ).jld")
state = ld["chain"][end]

# compute the coefficients that are going to be prescribed
H = create_FD_laplacian_reference_measure(q, μ).H
prescribed_values = H[:, end - 4^p + 1 : end]' * state

# defining the new sampler
cond_δ = 0.25
cond_k_max = k_max / 50



# creates sampler for density 
sampler = PCNMHSampler(create_FD_laplacian_reference_measure(q, μ, prescribed_values),
                       CahnHilliardDensity(ν),
                       cond_δ)


@time chain, acceptance_probability = sample_chain(sampler, cond_k_max, round(Int, cond_k_max / 10))

@show acceptance_probability

save("out/conditional_chain_p_$(p)_kmax_$(k_max)_nu_$(ν)_mu_$(μ)_T_$(T)_q_$(q)_delta_$(δ).jld", 
     "chain", chain,
     "μ", μ,
     "ν", ν,
     "T", T,
     "k_max", k_max,
     "q", q,
     "δ", δ)