using Agents

@agent struct CancerAgent(GridAgent{2})
    cell_cycle::Int
    age::Int
    dead::Bool
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
    if !agent.dead
        open_positions = collect(empty_nearby_positions(agent, model))

        if !isempty(open_positions)
            agent.age += 1
            if agent.age >= agent.cell_cycle
                if death(agent, model)
                    # agent dies
                    agent.dead = true
                    agent.cell_cycle = model.abtosis
                    agent.age = 0
                    if agent.age >= agent.cell_cycle
                        remove_agent!(agent, model)
                    end
                else
                    # agent proliferates
                    pos = rand(abmrng(model), open_positions)
                    
                    # 10% chance cell cycle mutates by ±1
                    new_cycle = agent.cell_cycle
                    r = rand(abmrng(model))
                    if r < model.evolution_rate / 2
                        new_cycle = max(model.t_min, agent.cell_cycle - 1)
                    elseif r < model.evolution_rate
                        new_cycle = min(model.t_max, agent.cell_cycle + 1)
                    end
                    
                    new_agent = CancerAgent(model, pos, new_cycle, 0, false)
                    add_agent_own_pos!(new_agent, model)
                end
            elseif rand(abmrng(model)) < model.velocity
                # agent moves
                pos = rand(abmrng(model), open_positions)
                move_agent!(agent, pos, model)
            end
        end
    else
        # dead agent countdown to removal
        agent.age += 1
        if agent.age >= agent.cell_cycle
            remove_agent!(agent, model)
        end
    end
end



function model_step!(model)
    if abmtime(model) % model.dosage_interval == 0
        model.treatment_function(model)
    end
end