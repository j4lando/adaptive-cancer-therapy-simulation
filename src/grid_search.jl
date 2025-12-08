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
@from "treatment.jl" import decreasing_adaptive_treatment!, smart_vacation_treatment!, smooth_adaptive_treatment!, default_adaptive_treatment!
@from "model.jl" import step!

function fluctuation_metric(N_history)
    return std(N_history)   # lower = more stable
end

function run_sim(alpha, beta, gamma;
                 steps=19000,
                 n0=3000, 
                 seeds=[42],
                 dosage_interval=10, di2=80,
                 max_cells=7000,
                 min_cells=2000,
                 c_mean = 1.0,
                 c_std = 0.25,
                 treatment_function = default_adaptive_treatment!)

    # Run multiple seeds and average results
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
            model.treatment_function = treatment_function
        )


        # Pre-allocate vector with known size
        N_hist = Vector{Int}(undef, steps)
        actual_steps = steps

        for t in 1:steps
            step!(model)
            n = nagents(model)
            N_hist[t] = n
            
            # Break if population exceeds bounds
            if n > max_cells || n < min_cells
                actual_steps = t
                break
            end
        end

        # Trim to actual steps taken
        N_hist = N_hist[1:actual_steps]
        
        fluctuations[s] = fluctuation_metric(N_hist)
        final_steps[s] = actual_steps
    end
    
    return (fluct=mean(fluctuations), survival_time=mean(final_steps)/24, fluct_std=std(fluctuations))
end

# Function to create heatmap for a specific gamma value
function create_heatmap(results, gamma_value; filename="heatmap_gamma_$(gamma_value).png")
    # Filter results for this gamma value
    filtered = filter(x -> x.gamma == gamma_value, results)
    
    if isempty(filtered)
        println("No results found for gamma = $gamma_value")
        return
    end
    
    # Get unique alpha and beta values
    unique_alphas = sort(unique([r.alpha for r in filtered]))
    unique_betas = sort(unique([r.beta for r in filtered]))
    
    # Create matrix for heatmap
    survival_matrix = zeros(length(unique_betas), length(unique_alphas))
    
    for r in filtered
        i = findfirst(==(r.beta), unique_betas)
        j = findfirst(==(r.alpha), unique_alphas)
        survival_matrix[i, j] = r.survival_time
    end
    
    # Create heatmap
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

alphas = 0:0.1:1.0
betas  = 0:0.1:1.0
gammas = 0:1:4
seeds = [42, 43, 44]  # Array of seeds to average over

params = [(α, β, γ) for α in alphas, β in betas, γ in gammas]
results = Vector{Any}(undef, length(params))

# Create atomic counter for thread-safe progress tracking
completed = Atomic{Int}(0)
total = length(params)

println("Starting grid search with $total parameter combinations (averaging over $(length(seeds)) seeds)...")

Threads.@threads for i in eachindex(params)
    α, β, γ = params[i]
    result = run_sim(α, β, γ; seeds=seeds, dosage_interval = 12, di2 =120, min_cells=10, max_cells=15000, treatment_function = smart_vacation_treatment!)
    results[i] = (alpha=α, beta=β, gamma=γ, fluct=result.fluct, survival_time=result.survival_time, fluct_std=result.fluct_std)
    
    # Thread-safe increment and print progress
    count = atomic_add!(completed, 1)
    if count % 10 == 0 || count == total  # Print every 10 or at end
        println("Progress: $count/$total ($(round(100*count/total, digits=1))%)")
    end
end

println("Grid search complete!")

# sorted = sort(filter(x -> x.survival_time == 1500, results), by = x -> x.survival_time)
sorted = sort(results, by = x -> -x.survival_time)
best = first(sorted)
println("Best parameters: ", best)
println("Survival time: ", best.survival_time, " ± ", best.fluct_std)

create_heatmap(results, 0; filename="heatmap_smart_vacation_treatment.png")
