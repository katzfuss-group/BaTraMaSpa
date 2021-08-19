using Base: Float64
using Zygote: gradient

# Abstract density and automatic differentiation 

# abstract type for a density 
abstract type AbstractLogDensity end

abstract type AbstractLogDensityGradient end

struct AutomaticDensityGradient{Td<:AbstractLogDensity}
    density::Td
end 

function ∇(density::AbstractLogDensity)
    return AutomaticDensityGradient(density)
end

function (grad::AutomaticDensityGradient)(x) 
    return first(gradient(x -> grad.density(x), x))
end

struct CahnHilliardDensity<:AbstractLogDensity
    ν::Float64
end

function (density::CahnHilliardDensity)(x)
    return density.ν * sum((1 .- x.^2).^2)
end