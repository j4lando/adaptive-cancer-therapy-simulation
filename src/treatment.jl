using Agents

"""
The adaptive dosage function as specified in the paper
"""
function default_adaptive_treatment!(model)
    n = nagents(model)
    n_0 = model.initial_agents
    if n < n_0 * 0.5
        model.dosage = 0
    elseif n > (1 + model.beta) * n_0
        model.dosage = min(1.0, (1 + model.alpha) * model.last_dosage)
        model.last_dosage = model.dosage
    elseif n <= (1 - model.beta) * n_0
        model.dosage = min(1.0, (1 - model.alpha) * model.last_dosage)
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
    D = model.last_dosage * S * (1.0 + model.alpha * tanh((x - 1.0) / model.beta))

    model.dosage = min(1.0, D)
end

"""
Josh's smooth dosage function
"""
function decreasing_adaptive_treatment!(model)
    if abmtime(model) % model.di2 == 0
        n = nagents(model)
        n_0 = model.initial_agents
        x = n / n_0

        # Smooth "on" function near 0.5*N0
        S = 1.0 / (1.0 + exp((0.5 - x) / model.beta))

        # Smooth modulation around N0
        D = S * (1.0 + model.alpha * tanh((x - 1.0) / model.beta))

        last_dose = max(0, 1 + min(0, 2*D - 2))
        first_dose = min(1, 2 - min(0, 2*D - 2))
        model.dosage = first_dose
        model.last_dosage = (last_dose - first_dose) / model.di2 * model.dosage_interval
    else
        model.dosage += model.last_dosage # last_dosage is an increment here
    end
end

