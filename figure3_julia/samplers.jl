# an abstract type for a sampler
abstract type AbstractSampler end 

# an abstract type for a basic MH sampler
struct PCNMHSampler{Tr,Td}
    reference_measure::Tr
    density::Td
    δ::Float64
    proposal::Vector{Float64}
end

# Constructor that automatically creates the scratch spaces
function PCNMHSampler(reference, density, δ)
    N = length(sample!(reference))
    return PCNMHSampler{typeof(reference), typeof(density)}(reference, 
                                                            density, 
                                                            δ, 
                                                            zeros(N))
end

function ρ(u, v, sampler::PCNMHSampler)
    return sampler.density(u)
end

function step!(state, sampler::PCNMHSampler)
    β = √(8 * sampler.δ) / (2 + sampler.δ)
    v = sampler.proposal
    u = state
    w = sample!(sampler.reference_measure)
    # The proposal v is obtained as v = √((1 - β^2)) * u + β * w,
    # where w is sampled according to the reference measure and u is the present state
    # proposal now contains a sample from the reference density 
    @. v = √((1 - β^2)) * u + β * w

    prescribed_H = sampler.reference_measure.H[:, (length(u) - size(sampler.reference_measure.A21, 2) + 1) : length(u)]
    v .+= prescribed_H * (prescribed_H' * (u - v))
     
    if exp(ρ(u, v, sampler) - ρ(v, u, sampler)) > rand()
        u .= v
        # returns true to signify that the proposal was accepted
        return true
    else 
        # returns false to signify that the proposal was rejected
        return false
    end
end