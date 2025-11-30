import Pkg; Pkg.add("Agents")
using Agents

@agent struct CancerAgent(GridAgent{2})
    cell_cycle::Int
    age::Int
end

function death(agent, model)
    sensitivity = 1 - (agent.cell_cycle - model.t_min) / (model.t_max - model.t_min)
    hill_coefficient = model.dosage ^ model.hill_n / (model.dosage ^ model.hill_n + model.hill_k ^ model.hill_n)

    if rand(abmrng(model)) < sensitivity * hill_coefficient
        return true
    end
    return false
end

function cancer_step!(agent, model)
    open_positions = collect(empty_nearby_positions(agent, model))

    if !isempty(open_positions)
        agent.age += 1
        if agent.age >= agent.cell_cycle
            if death(agent, model)
                # agent dies
                remove_agent!(agent, model)
            else
                # agent proliferates
                pos = rand(abmrng(model), open_positions)
                
                # 10% chance cell cycle mutates by ±1
                new_cycle = agent.cell_cycle
                r = rand(abmrng(model))
                if r < 0.05
                    new_cycle = max(model.t_min, agent.cell_cycle - 1)
                elseif r < 0.10
                    new_cycle = min(model.t_max, agent.cell_cycle + 1)
                end
                
                new_agent = CancerAgent(model, pos, new_cycle, 0)
                add_agent_own_pos!(new_agent, model)
            end
        elseif rand(abmrng(model)) < model.velocity
            # agent moves
            pos = rand(abmrng(model), open_positions)
            move_agent!(agent, pos, model)
        end
    end
end

function adjust_dosage!(model)
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

function model_step!(model)
    model.time += 1
    if model.time % model.dosage_interval == 0
        adjust_dosage!(model)
    end
    step!(model)
end