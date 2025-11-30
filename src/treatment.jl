import Pkg; Pkg.add("Agents")
using Agents

include("model.jl")

"""
The adaptive dosage function as specified in the paper
"""
function default_adaptive_treatment!(model)
    n = nagents(model)
    n_0 = model.initial_agents
    if n < n_0 * 0.5
        model.dosage = 0
    elseif n > (1 + model.beta) * n_0
        model.dosage = (1 + model.alpha) * model.last_dosage
        model.last_dosage = model.dosage
    elseif n <= (1 - model.beta) * n_0
        model.dosage = (1 - model.alpha) * model.last_dosage
        model.last_dosage = model.dosage
    end
end

"""
Li's smooth dosage function
"""
function smooth_adaptive_treatment!(model)
    n = nagents(model)
    n_0 = model.initial_agents
    x = n / n_0

    # Smooth "on" function near 0.5*N0
    S = 1.0 / (1.0 + exp((0.5 - x) / model.beta))

    # Smooth modulation around N0
    D = model.dosage * S * (1.0 + model.alpha * tanh((x - 1.0) / model.beta))

    model.dosage = D
end