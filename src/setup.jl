import Pkg
Pkg.add("Agents")
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("StatsBase")

using Agents
using Random
using Distributions
using StatsBase

include("model.jl")
include("treatment.jl")

function centered_positions(rng::AbstractRNG, cx::Int, cy::Int, σ::Float64, rmax::Float64, N::Int)
    # 1. Generate all integer grid points within rmax
    rmax_ceil = ceil(Int, rmax)
    candidates = Tuple{Int,Int}[]
    for x in (cx - rmax_ceil):(cx + rmax_ceil)
        for y in (cy - rmax_ceil):(cy + rmax_ceil)
            if sqrt((x - cx)^2 + (y - cy)^2) <= rmax
                push!(candidates, (x, y))
            end
        end
    end

    total_candidates = length(candidates)
    if N > total_candidates
        error("Cannot generate $N unique points: only $total_candidates available within radius $rmax")
    end

    # 2. Compute probabilities proportional to 2D normal distribution
    probs = [pdf(Normal(0, σ), x - cx) * pdf(Normal(0, σ), y - cy) for (x, y) in candidates]
    # Normalize to sum to 1
    probs ./= sum(probs)

    # 3. Sample N unique points according to probabilities
    chosen_indices = sample(rng, 1:total_candidates, Weights(probs), N; replace=false)
    return candidates[chosen_indices]
end

function initialize(
    num_agents, 
    seed; 
    dosage = 1.0, 
    hill_n = 1.5, 
    hill_k = 0.25, 
    t_min = 30, 
    t_max = 70, 
    size = 140, 
    move = false, 
    initial_variance = 10.0, 
    alpha = 0.25, 
    beta = 0.05, 
    velocity = 0.1,
    birth_population = 20,
    dosage_interval = 72,
    tail_skew = 0.3,
    treatment_function = default_adaptive_treatment!
)
    properties = Dict(
        :initial_agents => num_agents,
        :dosage => 0,
        :last_dosage => dosage,
        :hill_n => hill_n,
        :hill_k => hill_k,
        :t_min => t_min,
        :t_max => t_max,
        :move => move,
        :alpha => alpha,
        :beta => beta,
        :velocity => 0.05,
        :time => 0,
        :dosage_interval => dosage_interval,
        :treatment_function! => treatment_function
    )

    space = GridSpaceSingle((size, size), periodic = false, metric = :chebyshev)
    rng = Xoshiro(seed)

    model = StandardABM(
        CancerAgent, space; 
        agent_step! = cancer_step!, properties, rng)


    max_distance = sqrt(birth_population/3)
    center = Int(size / 2)
    initial_pos = centered_positions(rng, center, center, initial_variance, max_distance, birth_population)

    for pos in initial_pos
        distance = sqrt((pos[1] - center)^2 + (pos[2] - center)^2)
        normalized_distance = distance / max_distance
        
        # Apply tail skew: higher values push towards extremes (tmin/tmax)
        # tail_skew = 1.0 is linear, > 1.0 skews to tails, < 1.0 skews to center
        bias = 1 - normalized_distance^tail_skew
        
        cycle = round(Int, t_min + (t_max - t_min) * bias)
        
        add_agent!(pos, model; cell_cycle = cycle, age = 0)
        #add_agent!(pos, model; cell_cycle = rand(rng, t_min:t_max), age = 0)
    end

    condition(model, time) = nagents(model) >= num_agents
    step!(model, condition)
    model.velocity = velocity
    model.dosage = dosage
    model.initial_agents = nagents(model)

    return model
end
