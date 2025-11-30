import Pkg
Pkg.add("Agents")
Pkg.add("ColorSchemes")
Pkg.add("CairoMakie")

using Agents
using ColorSchemes
using CairoMakie

include("setup.jl")
include("model.jl")

function video(model; destination = "cancer.mp4", time_steps = 1000, break_condition = 13000)

    groupcolor(agent) = get(
        colorschemes[:rainbow], 
        (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    )

    fig, ax, abmobs = abmplot(
        model;
        agent_color = groupcolor, 
        agent_marker = :circle, 
        as = 10,
        title = "Cancer cell simulation"
    )

    record(fig, destination; framerate = 15) do io
        for i in 1:time_steps  # Maximum frames
            model_step!(model)
            abmobs.model[] = model  # Update the observable
            recordframe!(io)
            
            # Custom stopping condition
            if nagents(model) > break_condition
                break
            end
        end
    end
end