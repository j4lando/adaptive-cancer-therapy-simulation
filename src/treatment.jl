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
        x = n / model.gamma / n_0

        # # Smooth "on" function near 0.5*N0
        # S = 1.0 / (1.0 + exp((0.5 - x) / model.beta))

        # # Smooth modulation around N0
        # D = S * (1.0 + model.alpha * tanh((x - 1.0) / model.beta))
        D = model.alpha + model.beta * x

        last_dose = max(0, 1 + min(0, 2*D - 2))
        first_dose = min(1, 2 + min(0, 2*D - 2))
        model.dosage = first_dose
        model.last_dosage = (last_dose - first_dose) / model.di2 * model.dosage_interval
    else
        model.dosage = clamp(model.dosage + model.last_dosage, 0.0, 1.0)  # ← Clamp to [0,1]
    end
end

function smart_vacation_treatment!(model)
    if abmtime(model) % model.di2 == 0
        n = nagents(model)
        n_0 = model.initial_agents
        dosage = model.dosage + model.last_dosage
        if n > (1 + model.beta) * n_0
            dosage = min(2, (1 + model.alpha) * dosage)
        elseif n <= (1 - model.beta) * n_0
            dosage = min(2, (1 - model.alpha) * dosage)
        end
        model.dosage = clamp(dosage, 0.0, 1.0)
        model.last_dosage = clamp((dosage - 1.0), 0.0, 1.0)
    else
        model.dosage, model.last_dosage = model.last_dosage, model.dosage
    end
end

