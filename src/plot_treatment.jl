
using CairoMakie

"""
Compute the default adaptive treatment dosage based on population n relative to baseline n_0.
Returns the new dosage and last_dosage values.
"""
function default_adaptive_treatment(n, n_0, alpha, beta, initial_dosage = 1.0)  
    if n < n_0 * 0.5
        new_dosage = 0
    elseif n > (1 + beta) * n_0
        new_dosage = (1 + alpha) * initial_dosage
    elseif n <= (1 - beta) * n_0
        new_dosage = (1 - alpha) * initial_dosage
    else
        new_dosage = initial_dosage
    end

    return min(new_dosage, 1.0)
end

"""
Compute the smooth adaptive treatment dosage using continuous functions.
"""
function smooth_adaptive_treatment(n, n_0, alpha, beta, initial_dosage = 1.0)
    x = n / n_0
    
    # Smooth "on" function near 0.5*N0
    S = 1.0 / (1.0 + exp((0.5 - x) / beta))
    
    # Smooth modulation around N0
    D = initial_dosage * S * (1.0 + alpha * tanh((x - 1.0) / beta))
    
    return min(D, 1.0)
end

"""
Visualize how treatment dosage changes with population for both strategies.
"""
function plot_treatment_strategies(; 
    n_0 = 3000, 
    alpha = 0.25, 
    beta = 0.05, 
    initial_dosage = 1.0,
    n_max = 6000,
    filename = "treatment_comparison.png"
)
    # Create a range of population values
    n_range = 0:10:n_max
    
    # Calculate dosages for default adaptive treatment
    default_dosages = [default_adaptive_treatment(n, n_0, alpha, beta, initial_dosage) for n in n_range]
    
    # Calculate dosages for smooth adaptive treatment
    smooth_dosages = [smooth_adaptive_treatment(n, n_0, alpha, beta, initial_dosage) for n in n_range]
    
    # Create the plot
    fig = Figure(size = (1000, 600))
    
    ax = Axis(fig[1, 1],
        xlabel = "Population (n)",
        ylabel = "Dosage",
        title = "Adaptive Treatment Strategies Comparison",
        xlabelsize = 14,
        ylabelsize = 14,
        titlesize = 16
    )
    
    # Plot both treatment strategies
    lines!(ax, collect(n_range), default_dosages, 
        label = "Default Adaptive", linewidth = 2.5, color = :blue)
    lines!(ax, collect(n_range), smooth_dosages, 
        label = "Smooth Adaptive", linewidth = 2.5, color = :red)
    
    # Add reference lines for key thresholds
    vlines!(ax, [n_0 * 0.5], 
        color = :gray, linestyle = :dash, linewidth = 1.5, label = "0.5×N₀ (shutdown)")
    vlines!(ax, [n_0], 
        color = :black, linestyle = :dot, linewidth = 1.5, label = "N₀ (baseline)")
    vlines!(ax, [(1 + beta) * n_0], 
        color = :orange, linestyle = :dash, linewidth = 1.5, label = "(1+β)×N₀ (increase)")
    vlines!(ax, [(1 - beta) * n_0], 
        color = :green, linestyle = :dash, linewidth = 1.5, label = "(1-β)×N₀ (decrease)")
    
    axislegend(ax, position = :lt)
    
    # Save and display
    save(filename, fig)
    println("Plot saved as $filename")
    
    return fig
end

# Run the visualization
plot_treatment_strategies(initial_dosage = 1.0)
