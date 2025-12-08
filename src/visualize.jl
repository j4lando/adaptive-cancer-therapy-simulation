using Agents
using ColorSchemes
using CairoMakie
using Printf
using DataFrames

import FromFile: @from
@from "model.jl" import model_step! 

function video(model; destination = "cancer.mp4", time_steps = 1000, break_condition = 13000, fps = 30, resolution = (1000, 800), time_resolution = 1)

    # White (0% resistance) to Red (100% resistance)
    resistance_cmap = cgrad(:Reds, rev = false)
    groupcolor(agent) = get(
        resistance_cmap, 
        (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    )

    fig = Figure(resolution = resolution)
    
    # Main ABM plot area
    ax_abm = Axis(fig[1, 1], title = "Cancer cell simulation")
    
    # Create colorbar for resistance
    Colorbar(fig[1, 2], 
        colormap = resistance_cmap,
        limits = (0, 1),
        label = "Resistance",
        width = 20)
    
    # Dosage indicator text
    dosage_label = Label(fig[2, 1:2], 
        text = @sprintf("Current Dosage: %.3f", model.dosage),
        tellwidth = false,
        fontsize = 16,
        halign = :center)

    abmobs = abmplot!(
        ax_abm,
        model;
        agent_color = groupcolor, 
        agent_marker = :circle, 
        as = 10
    )
    
    # Initialize vectors for tracking metrics (much faster than DataFrame)
    dosages = Vector{Float64}(undef, time_steps + 1)
    agent_counts = Vector{Int}(undef, time_steps + 1)
    
    # Record initial state
    dosages[1] = model.dosage
    agent_counts[1] = nagents(model)
    
    actual_steps = 1

    record(fig, destination; framerate = fps) do io
        for i in 1:time_steps
            step!(model)
            abmobs.model[] = model
            
            # Update metrics tracking
            dosages[i + 1] = model.dosage
            agent_counts[i + 1] = nagents(model)
            
            # Update dosage label
            dosage_label.text = @sprintf("Current Dosage: %.3f | Agents: %d", model.dosage, nagents(model))
            
            if i % time_resolution == 0
                recordframe!(io)
            end
            
            
            actual_steps = i + 1
            
            # Custom stopping condition
            if nagents(model) > break_condition || nagents(model) == 0
                break
            end
        end
    end
    
    # Trim vectors to actual steps taken
    resize!(dosages, actual_steps)
    resize!(agent_counts, actual_steps)
    
    # Create time linespace starting at 0
    time_points = range(0, actual_steps - 1, length = actual_steps)
    
    # Save dosage plot as separate graph
    fig_dosage = Figure(resolution = (800, 600))
    ax_dosage_only = Axis(fig_dosage[1, 1],
        title = "Dosage over Time",
        xlabel = "Time",
        ylabel = "Dosage")
    lines!(ax_dosage_only, time_points, dosages, color = :blue, linewidth = 2)
    ylims!(ax_dosage_only, 0, 1.5)
    
    dosage_graph_path = replace(destination, ".mp4" => "_dosage.png")
    save(dosage_graph_path, fig_dosage)
    
    # Save agent count plot as separate graph
    fig_agents = Figure(resolution = (800, 600))
    ax_agents = Axis(fig_agents[1, 1],
        title = "Agent Count over Time",
        xlabel = "Time",
        ylabel = "Number of Agents")
    lines!(ax_agents, time_points, agent_counts, color = :red, linewidth = 2)
    
    agents_graph_path = replace(destination, ".mp4" => "_agents.png")
    save(agents_graph_path, fig_agents)
end


function snapshot(model; destination = "cancer_snapshot.png", resolution = (1000, 800))
    
    # White (0% resistance) to Red (100% resistance)
    resistance_cmap = cgrad(:Reds, rev = false)
    groupcolor(agent) = get(
        resistance_cmap, 
        (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    )

    fig = Figure(resolution = resolution)
    
    # Main plot area
    ax = Axis(fig[1, 1], title = "Cancer cell simulation")
    
    # Create colorbar for resistance
    Colorbar(fig[1, 2], 
        colormap = resistance_cmap,
        limits = (0, 1),
        label = "Resistance",
        width = 20)
    
    # Dosage and agent count info
    Label(fig[2, 1:2], 
        text = @sprintf("Timestep: %d | Dosage: %.3f | Agents: %d", 
            abmtime(model), model.dosage, nagents(model)),
        tellwidth = false,
        fontsize = 16,
        halign = :center)

    abmplot!(
        ax,
        model;
        agent_color = groupcolor, 
        agent_marker = :circle, 
        as = 10
    )
    
    save(destination, fig)
    return fig
end

function multi_snapshot(models; destinations = nothing, time_steps = 1000, break_condition = 13000, interval = 100, combined_name = "combined")
    
    # Set default destinations if not provided
    if destinations === nothing
        destinations = ["model_$(i)" for i in 1:length(models)]
    end
    
    @assert length(models) == length(destinations) "Number of models must match number of destinations"
    
    num_models = length(models)
    
    # Initialize storage for all models
    all_dosages = [Vector{Float64}(undef, time_steps + 1) for _ in 1:num_models]
    all_agent_counts = [Vector{Int}(undef, time_steps + 1) for _ in 1:num_models]
    all_actual_steps = ones(Int, num_models)
    snapshot_counts = zeros(Int, num_models)
    
    # Record initial states and save initial snapshots
    for (m_idx, (model, destination)) in enumerate(zip(models, destinations))
        # White (0% resistance) to Red (100% resistance)
        resistance_cmap = cgrad(:Reds, rev = false)
        groupcolor(agent) = get(
            resistance_cmap, 
            (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
        )
        
        all_dosages[m_idx][1] = model.dosage
        all_agent_counts[m_idx][1] = nagents(model)
        
        # Save initial snapshot
        fig = Figure(resolution = (1000, 800))
        ax = Axis(fig[1, 1], title = "Cancer cell simulation - Step 0")
        
        Colorbar(fig[1, 2], 
            colormap = resistance_cmap,
            limits = (0, 1),
            label = "Resistance",
            width = 20)
        
        Label(fig[2, 1:2], 
            text = @sprintf("Timestep: 0 | Dosage: %.3f | Agents: %d", 
                model.dosage, nagents(model)),
            tellwidth = false,
            fontsize = 16,
            halign = :center)
        
        abmplot!(ax, model; agent_color = groupcolor, agent_marker = :circle, as = 10)
        save("$(destination)_$(lpad(snapshot_counts[m_idx], 4, '0')).png", fig)
        snapshot_counts[m_idx] += 1
    end
    
    # Run simulation for all models
    for i in 1:time_steps
        for (m_idx, (model, destination)) in enumerate(zip(models, destinations))
            # Skip this model if it already stopped
            if nagents(model) > break_condition || nagents(model) == 0
                continue
            end
            
            step!(model)
            
            # Update metrics tracking
            all_dosages[m_idx][i + 1] = model.dosage
            all_agent_counts[m_idx][i + 1] = nagents(model)
            
            all_actual_steps[m_idx] = i + 1
            
            # Take snapshot at intervals
            if i % interval == 0
                # White (0% resistance) to Red (100% resistance)
                resistance_cmap = cgrad(:Reds, rev = false)
                groupcolor(agent) = get(
                    resistance_cmap, 
                    (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
                )
                
                fig = Figure(resolution = (1000, 800))
                ax = Axis(fig[1, 1], title = "Cancer cell simulation - Step $i")
                
                Colorbar(fig[1, 2], 
                    colormap = resistance_cmap,
                    limits = (0, 1),
                    label = "Resistance",
                    width = 20)
                
                Label(fig[2, 1:2], 
                    text = @sprintf("Timestep: %d | Dosage: %.3f | Agents: %d", 
                        i, model.dosage, nagents(model)),
                    tellwidth = false,
                    fontsize = 16,
                    halign = :center)
                
                abmplot!(ax, model; agent_color = groupcolor, agent_marker = :circle, as = 10)
                save("$(destination)_$(lpad(snapshot_counts[m_idx], 4, '0')).png", fig)
                snapshot_counts[m_idx] += 1
            end
        end
        
        # Break if all models have stopped
        if all(m_idx -> nagents(models[m_idx]) > break_condition || nagents(models[m_idx]) == 0, 1:num_models)
            break
        end
    end
    
    # Process each model's data and create individual dosage plots
    for (m_idx, destination) in enumerate(destinations)
        actual_steps = all_actual_steps[m_idx]
        
        # Trim vectors to actual steps taken
        dosages = all_dosages[m_idx][1:actual_steps]
        
        # Convert time steps to days
        time_days = range(0, (actual_steps - 1) / 24, length = actual_steps)
        
        # Save dosage plot with filled area
        fig_dosage = Figure(resolution = (800, 600))
        ax_dosage = Axis(fig_dosage[1, 1],
            title = "Dosage over Time",
            xlabel = "Time (days)",
            ylabel = "Dosage")
        
        # Fill area under curve
        band!(ax_dosage, time_days, zeros(length(dosages)), dosages, color = (:blue, 0.3))
        lines!(ax_dosage, time_days, dosages, color = :blue, linewidth = 2)
        ylims!(ax_dosage, 0, 1.5)
        save("$(destination)_dosage.png", fig_dosage)
        
        println("Model $(m_idx): Saved $(snapshot_counts[m_idx]) snapshots and dosage graph")
    end
    
    # Create combined agent count plot
    fig_agents = Figure(resolution = (800, 600))
    ax_agents = Axis(fig_agents[1, 1],
        title = "Agent Count over Time (All Models)",
        xlabel = "Time (days)",
        ylabel = "Number of Agents")
    
    colors = [:red, :green, :purple, :orange, :cyan, :magenta]
    
    for (m_idx, destination) in enumerate(destinations)
        actual_steps = all_actual_steps[m_idx]
        agent_counts = all_agent_counts[m_idx][1:actual_steps]
        time_days = range(0, (actual_steps - 1) / 24, length = actual_steps)
        
        color = colors[mod1(m_idx, length(colors))]
        lines!(ax_agents, time_days, agent_counts, color = color, linewidth = 4, label = destination)
    end
    
    ylims!(ax_agents, 0, break_condition)
    axislegend(ax_agents, position = :lt)
    save("$(combined_name)_agents.png", fig_agents)
    
    println("Saved combined agent count graph as $(combined_name)_agents.png")
end

