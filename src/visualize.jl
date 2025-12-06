using Agents
using ColorSchemes
using CairoMakie
using Printf

import FromFile: @from
@from "model.jl" import model_step! 

function video(model; destination = "cancer.mp4", time_steps = 1000, break_condition = 13000)

    groupcolor(agent) = get(
        colorschemes[:rainbow], 
        (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    )

    fig = Figure(resolution = (1000, 800))
    
    # Main plot area
    ax = Axis(fig[1, 1], title = "Cancer cell simulation")
    
    # Create colorbar for cell cycle
    colorbar_ax = Axis(fig[1, 2], 
        width = 20,
        ylabel = "Cell Cycle Time",
        ylabelrotation = 0,
        ylabelsize = 12)
    
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
        ax,
        model;
        agent_color = groupcolor, 
        agent_marker = :circle, 
        as = 10
    )

    record(fig, destination; framerate = 70) do io
        for i in 1:time_steps  # Maximum frames
            model_step!(model)
            abmobs.model[] = model  # Update the observable
            
            # Update dosage label
            dosage_label.text = @sprintf("Current Dosage: %.3f | Agents: %d", model.dosage, nagents(model))
            
            recordframe!(io)
            
            # Custom stopping condition
            if nagents(model) > break_condition
                break
            end
        end
    end
end