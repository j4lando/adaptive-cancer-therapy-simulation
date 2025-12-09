using Agents
using Random: Xoshiro
using Distributions
using StatsBase
using ColorSchemes
using CairoMakie
using Printf
using Statistics: mean

import FromFile: @from 
@from "setup.jl" import initialize
@from "visualize.jl" import video, snapshot, multi_snapshot
@from "treatment.jl" import default_adaptive_treatment!, smooth_adaptive_treatment!, smart_vacation_treatment! 

# ============================================================================
# Simulation Constants
# ============================================================================

const EVOLUTION_RATE = 0.1
const VELOCITY = 0.1
const INITIAL_VELOCITY = 0.01
const BIRTH_POPULATION = 100
const GRID_SIZE = 200
const INITIAL_AGENTS = 15000
const ABTOSIS = 20
const DOSAGE = 0.5

# ============================================================================
# Experiment 1: Treatment Comparison (New)
# ============================================================================

c_means = [1.0, 1.0, 0.6, 0.6]
c_stds = [0.05, 0.25, 0.05, 0.25]
prefixes = ["A_new", "B_new", "C_new", "D_new"]
seed = 70

for i in 1:length(c_means)
    c_mean = c_means[i]
    c_std = c_stds[i]
    prefix = prefixes[i]
    
    model_smooth = initialize(INITIAL_AGENTS, seed+i, alpha=4.4, beta=1.2, 
        birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=DOSAGE, 
        treatment_function=smooth_adaptive_treatment!)
    
    modelat2 = initialize(INITIAL_AGENTS, seed+i, alpha=0.5, beta=0.1, 
        birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=DOSAGE)
    
    model_slow_vacation = initialize(INITIAL_AGENTS, seed+i, dosage_interval=12, di2=120, 
        alpha=0.1, beta=0.8, birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=DOSAGE, 
        treatment_function=smart_vacation_treatment!)
    
    model_fast_vacation = initialize(INITIAL_AGENTS, seed+i, dosage_interval=12, di2=72, 
        alpha=0.5, beta=0.5, birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=DOSAGE, 
        treatment_function=smart_vacation_treatment!)
    
    modelat3 = initialize(INITIAL_AGENTS, seed+i, alpha=0.9, beta=1.0, 
        birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=DOSAGE)
    


    multi_snapshot(
        [model_smooth, modelat2, model_slow_vacation, model_fast_vacation, modelat3];
        destinations = ["$(prefix)_treatment_SMO", "$(prefix)_treatment_AT2", "$(prefix)_treatment_SV", "$(prefix)_treatment_FV", "$(prefix)_treatment_AT3"],
        time_steps = 19000,
        interval = 1000,
        break_condition = 35000,
        combined_name = "$(prefix)_treatment_comparison"
    )
end

# ============================================================================
# Experiment 2: Treatment Comparison (Original)
# ============================================================================

c_means = [1.0, 1.0, 0.6, 0.6]
c_stds = [0.05, 0.25, 0.05, 0.25]
prefixes = ["A", "B", "C", "D"]

for i in 1:length(c_means)
    c_mean = c_means[i]
    c_std = c_stds[i]
    prefix = prefixes[i]
    
    modelat1 = initialize(INITIAL_AGENTS, 20, 
        birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=DOSAGE, abtosis=ABTOSIS)
    
    modelat2 = initialize(INITIAL_AGENTS, 20, alpha=0.5, beta=0.1, 
        birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=DOSAGE, abtosis=ABTOSIS)
    
    modelct = initialize(INITIAL_AGENTS, 20, dosage_interval=1000000, 
        birth_population=BIRTH_POPULATION, mean=c_mean, std=c_std, 
        evolution_rate=EVOLUTION_RATE, initial_velocity=INITIAL_VELOCITY, 
        velocity=VELOCITY, size=GRID_SIZE, dosage=1.0, abtosis=ABTOSIS)



    multi_snapshot(
        [modelat1, modelat2, modelct];
        destinations = ["$(prefix)_treatment_AT1", "$(prefix)_treatment_AT2", "$(prefix)_treatment_CT"],
        time_steps = 19000,
        interval = 4000,
        break_condition = 35000,
        combined_name = "$(prefix)_comparison"
    )
end
