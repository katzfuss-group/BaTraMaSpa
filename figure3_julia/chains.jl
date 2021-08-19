# samples a chain from a sampler over k steps
function sample_chain(sampler, k_max, warmup=10000)
    # creates the state from the reference measure
    state = copy(vec(sample!(sampler.reference_measure)))
    for k = 1 : warmup
        step!(state, sampler)
    end

    chain = typeof(state)[]
    push!(chain, copy(state))
    for k = 1 : k_max
        if step!(state, sampler)
            push!(chain, copy(state))
        end
        if mod(k, round(Int, k_max / 10)) == 0
            println("k = $(k)")
        end
    end
    return chain, (length(chain) - 1) / k_max
end

function sample_chain(sampler, k_max, state::AbstractVector, warmup=0)
    # creates the state from the reference measure
    for k = 1 : warmup
        step!(state, sampler)
    end

    chain = typeof(state)[]
    push!(chain, copy(state))
    for k = 1 : k_max
        if step!(state, sampler)
            push!(chain, copy(state))
        end
        if mod(k, round(Int, k_max / 10)) == 0
            println("k = $(k)")
        end
    end
    return chain, (length(chain) - 1) / k_max
end

