using Agents
using ColorSchemes
using CairoMakie
using Printf
using DataFrames

import FromFile: @from
@from "model.jl" import model_step! 

function video(model; destination = "cancer.mp4", time_steps = 1000, break_condition = 13000, fps = 30, resolution = (1000, 800))

    groupcolor(agent) = get(
        colorschemes[:rainbow], 
        (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    )

    fig = Figure(resolution = resolution)
    
    # Main ABM plot area
    ax_abm = Axis(fig[1, 1], title = "Cancer cell simulation")
    
    # Create colorbar for cell cycle
    Colorbar(fig[1, 2], 
        colormap = :rainbow,
        limits = (model.t_min, model.t_max),
        label = "Cell Cycle Time",
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
            
            recordframe!(io)
            
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
    
    groupcolor(agent) = get(
        colorschemes[:rainbow], 
        (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    )

    fig = Figure(resolution = resolution)
    
    # Main plot area
    ax = Axis(fig[1, 1], title = "Cancer cell simulation")
    
    # Create colorbar for cell cycle
    Colorbar(fig[1, 2], 
        colormap = :rainbow,
        limits = (model.t_min, model.t_max),
        label = "Cell Cycle Time",
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

function multi_snapshot(model; destination = "cancer_snapshots", time_steps = 1000, break_condition = 13000, interval = 100)
    
    groupcolor(agent) = get(
        colorschemes[:rainbow], 
        (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    )
    
    # Initialize vectors for tracking metrics
    dosages = Vector{Float64}(undef, time_steps + 1)
    agent_counts = Vector{Int}(undef, time_steps + 1)
    
    # Record initial state
    dosages[1] = model.dosage
    agent_counts[1] = nagents(model)
    
    # Save initial snapshot
    snapshot_count = 0
    fig = Figure(resolution = (1000, 800))
    ax = Axis(fig[1, 1], title = "Cancer cell simulation - Step 0")
    
    Colorbar(fig[1, 2], 
        colormap = :rainbow,
        limits = (model.t_min, model.t_max),
        label = "Cell Cycle Time",
        width = 20)
    
    Label(fig[2, 1:2], 
        text = @sprintf("Timestep: 0 | Dosage: %.3f | Agents: %d", 
            model.dosage, nagents(model)),
        tellwidth = false,
        fontsize = 16,
        halign = :center)
    
    abmplot!(ax, model; agent_color = groupcolor, agent_marker = :circle, as = 10)
    save("$(destination)_$(lpad(snapshot_count, 4, '0')).png", fig)
    snapshot_count += 1
    
    actual_steps = 1
    
    # Run simulation
    for i in 1:time_steps
        step!(model)
        
        # Update metrics tracking
        dosages[i + 1] = model.dosage
        agent_counts[i + 1] = nagents(model)
        
        actual_steps = i + 1
        
        # Take snapshot at intervals
        if i % interval == 0
            fig = Figure(resolution = (1000, 800))
            ax = Axis(fig[1, 1], title = "Cancer cell simulation - Step $i")
            
            Colorbar(fig[1, 2], 
                colormap = :rainbow,
                limits = (model.t_min, model.t_max),
                label = "Cell Cycle Time",
                width = 20)
            
            Label(fig[2, 1:2], 
                text = @sprintf("Timestep: %d | Dosage: %.3f | Agents: %d", 
                    i, model.dosage, nagents(model)),
                tellwidth = false,
                fontsize = 16,
                halign = :center)
            
            abmplot!(ax, model; agent_color = groupcolor, agent_marker = :circle, as = 10)
            save("$(destination)_$(lpad(snapshot_count, 4, '0')).png", fig)
            snapshot_count += 1
        end
        
        # Custom stopping condition
        if nagents(model) > break_condition || nagents(model) == 0
            break
        end
    end
    
    # Trim vectors to actual steps taken
    resize!(dosages, actual_steps)
    resize!(agent_counts, actual_steps)
    
    # Create time linespace starting at 0
    time_points = range(0, actual_steps - 1, length = actual_steps)
    
    # Save dosage plot
    fig_dosage = Figure(resolution = (800, 600))
    ax_dosage = Axis(fig_dosage[1, 1],
        title = "Dosage over Time",
        xlabel = "Time",
        ylabel = "Dosage")
    lines!(ax_dosage, time_points, dosages, color = :blue, linewidth = 2)
    ylims!(ax_dosage, 0, 1.5)
    save("$(destination)_dosage.png", fig_dosage)
    
    # Save agent count plot
    fig_agents = Figure(resolution = (800, 600))
    ax_agents = Axis(fig_agents[1, 1],
        title = "Agent Count over Time",
        xlabel = "Time",
        ylabel = "Number of Agents")
    lines!(ax_agents, time_points, agent_counts, color = :red, linewidth = 2)
    save("$(destination)_agents.png", fig_agents)
    
    println("Saved $snapshot_count snapshots and 2 metric graphs")
end

