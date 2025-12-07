using Agents
using Random
using Distributions
using StatsBase
using ColorSchemes
using CairoMakie
using Printf
using Statistics: mean, std
using Base.Threads: Atomic, atomic_add!

import FromFile: @from 
@from "setup.jl" import initialize
@from "treatment.jl" import decreasing_adaptive_treatment!
@from "model.jl" import step!

function fluctuation_metric(N_history)
    return std(N_history)   # lower = more stable
end

function run_sim(alpha, beta, gamma;
                 steps=1500,
                 n0=3000, seed=42,
                 dosage_interval=10, di2=80,
                 max_cells=13000)

    model = initialize(n0, seed;
        dosage_interval = dosage_interval,
        di2 = di2,
        dosage = 0.5,
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        birth_population = 100, tail_skew = .5, initial_velocity = .1, initial_variance = 3.0, velocity = 0.1
    )
    model.treatment_function = decreasing_adaptive_treatment!

    # Pre-allocate vector with known size
    N_hist = Vector{Int}(undef, steps)
    actual_steps = steps

    for t in 1:steps
        step!(model)
        n = nagents(model)
        N_hist[t] = n
        
        # Break if population exceeds max_cells
        if n > max_cells
            actual_steps = t
            break
        end
    end

    # Trim to actual steps taken
    N_hist = N_hist[1:actual_steps]
    
    return (fluct=fluctuation_metric(N_hist), final_cells=N_hist[end])
end

alphas = 0:0.1:1.0
betas  = 0:0.1:1.0
gammas = 0:1:5

params = [(α, β, γ) for α in alphas, β in betas, γ in gammas]
results = Vector{Any}(undef, length(params))

# Create atomic counter for thread-safe progress tracking
completed = Atomic{Int}(0)
total = length(params)

println("Starting grid search with $total parameter combinations...")

Threads.@threads for i in eachindex(params)
    α, β, γ = params[i]
    result = run_sim(α, β, γ)
    results[i] = (alpha=α, beta=β, gamma=γ, fluct=result.fluct, final_cells=result.final_cells)
    
    # Thread-safe increment and print progress
    count = atomic_add!(completed, 1)
    if count % 10 == 0 || count == total  # Print every 10 or at end
        println("Progress: $count/$total ($(round(100*count/total, digits=1))%)")
    end
end

println("Grid search complete!")

sorted = sort(filter(x -> x.final_cells > 2000, results), by = x -> x.fluct)
best = first(sorted)
println("Best parameters: ", best)
println("Final cell count: ", best.final_cells)
