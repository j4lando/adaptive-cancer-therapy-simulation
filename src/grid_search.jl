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
@from "treatment.jl" import decreasing_adaptive_treatment!, smart_vacation_treatment!, smooth_adaptive_treatment!
@from "model.jl" import step!

# ============================================================================
# Metrics and Simulation Functions
# ============================================================================

function fluctuation_metric(N_history)
    return std(N_history)
end

function run_sim(alpha, beta, gamma;
                 steps=19000,
                 n0=3000, 
                 seeds=[42],
                 dosage_interval=10, 
                 di2=80,
                 max_cells=7000,
                 min_cells=2000,
                 max_dosage_steps=4,
                 c_mean=1.0,
                 c_std=0.25,
                 size=150,
                 abtosis=20,
                 treatment_function=default_adaptive_treatment!)

    num_seeds = length(seeds)
    fluctuations = Vector{Float64}(undef, num_seeds)
    final_steps = Vector{Int}(undef, num_seeds)
    
    for (s, seed) in enumerate(seeds)
        model = initialize(n0, seed;
            dosage_interval = dosage_interval,
            di2 = di2,
            dosage = 0.5,
            alpha = alpha,
            beta = beta,
            gamma = gamma,
            birth_population = 100, 
            mean = c_mean, std = c_std, 
            evolution_rate = 0.1, 
            initial_velocity = .01, 
            velocity = 0.1, 
            size = 150,
            abtosis = abtosis,
            treatment_function = treatment_function
        )

        N_hist = Vector{Int}(undef, steps)
        actual_steps = steps
        max_dosage_counter = 0
        prev_n = n0

        for t in 1:steps
            step!(model)
            n = nagents(model)
            N_hist[t] = n
            
            # Check if both dosages are at maximum and population is increasing
            if model.dosage >= 1.0 && model.last_dosage >= 1.0 && n > prev_n
                max_dosage_counter += 1
            else
                max_dosage_counter = 0
            end
            
            # Break if population exceeds bounds or max dosage fails for too long
            if n > max_cells || n <= min_cells || max_dosage_counter >= max_dosage_steps * dosage_interval || n == size^2
                actual_steps = t
                break
            end
            
            prev_n = n
        end

        N_hist = N_hist[1:actual_steps]
        fluctuations[s] = fluctuation_metric(N_hist)
        final_steps[s] = actual_steps
    end
    
    return (fluct=mean(fluctuations), survival_time=mean(final_steps)/24, fluct_std=std(fluctuations))
end

# ============================================================================
# AI Generated Visualization Functions
# ============================================================================

function create_heatmap(results, gamma_value; filename="heatmap_gamma_$(gamma_value).png", clamp_min=nothing)
    # Filter results for this gamma value
    filtered = filter(x -> x.gamma == gamma_value, results)
    
    if isempty(filtered)
        println("No results found for gamma = $gamma_value")
        return
    end
    
    unique_alphas = sort(unique([r.alpha for r in filtered]))
    unique_betas = sort(unique([r.beta for r in filtered]))
    
    survival_matrix = zeros(length(unique_betas), length(unique_alphas))
    
    for r in filtered
        i = findfirst(==(r.beta), unique_betas)
        j = findfirst(==(r.alpha), unique_alphas)
        survival_matrix[i, j] = r.survival_time
    end
    
    if !isnothing(clamp_min)
        survival_matrix = max.(survival_matrix, clamp_min)
    end
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],
        xlabel = "Alpha",
        ylabel = "Beta",
        title = "Survival Time Heatmap")
    
    hm = heatmap!(ax, unique_alphas, unique_betas, survival_matrix',
        colormap = :viridis)
    
    Colorbar(fig[1, 2], hm, label = "Survival Time (steps)")
    
    save(filename, fig)
    println("Saved heatmap to $filename")
    
    return fig
end

# ============================================================================
# Grid Search Configuration and Execution
# ============================================================================

# Parameter ranges
alphas = 1.0:0.2:5.0
betas = 0.0:0.1:1.2
gammas = 0:1:0
seeds = [52, 53, 54]

# Simulation settings
SIM_CONFIG = (
    n0 = 15000,
    dosage_interval = 72,
    di2 = 72,
    min_cells = 0,
    max_cells = 35000,
    treatment_function = smooth_adaptive_treatment!,
    size = 200,
    abtosis = 20
)

# Generate all parameter combinations
params = [(α, β, γ) for α in alphas, β in betas, γ in gammas]
results = Vector{Any}(undef, length(params))

# Thread-safe progress tracking
completed = Atomic{Int}(0)
total = length(params)

println("Starting grid search with $total parameter combinations (averaging over $(length(seeds)) seeds)...")

Threads.@threads for i in eachindex(params)
    α, β, γ = params[i]
    result = run_sim(α, β, γ; seeds=seeds, SIM_CONFIG...)
    results[i] = (alpha=α, beta=β, gamma=γ, fluct=result.fluct, survival_time=result.survival_time, fluct_std=result.fluct_std)
    
    count = atomic_add!(completed, 1)
    if count % 10 == 0 || count == total
        println("Progress: $count/$total ($(round(100*count/total, digits=1))%)")
    end
end

println("Grid search complete!")

# ============================================================================
# Results Analysis
# ============================================================================

sorted = sort(results, by = x -> -x.survival_time)
best = first(sorted)
println("Best parameters: ", best)
println("Survival time: ", best.survival_time, " ± ", best.fluct_std)

name = "smooth_medium"
create_heatmap(results, 0; filename="heatmap_$(name)_treatment.png")
create_heatmap(results, 0; filename="heatmap_$(name)_treatment_higher_resolution.png", clamp_min=90)