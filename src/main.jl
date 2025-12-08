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
@from "treatment.jl" import default_adaptive_treatment!, smooth_adaptive_treatment!, decreasing_adaptive_treatment!, smart_vacation_treatment! 

cancer = initialize(3000, 42)
# snapshot(cancer, destination = "initial_snapshot.png")
# # video(cancer, time_steps = 1500, destination = "adaptive_treatment.mp4")

cancer = initialize(3000, 42)
multi_snapshot(cancer, time_steps = 1500, interval = 300, destination = "adaptive_treatment_multi")


cancer = initialize(3000, 44, dosage_interval = 72, di2 = 80, dosage = 0.5, alpha = 0.5, beta = .1, gamma = 1.0, birth_population = 100, tail_skew = .5, initial_velocity = 0.1, initial_variance = 3.0, velocity = 0.1, evolution_rate = 0.0)
multi_snapshot(cancer, time_steps = 9000, interval = 1000, destination = "paper_dosage")

cancer = initialize(3000, 44, dosage_interval = 72, di2 = 80, dosage = 0.5, alpha = 0.85, beta = .15, gamma = 1.0, birth_population = 100, tail_skew = .5, initial_velocity = 0.1, initial_variance = 3.0, velocity = 0.1, evolution_rate = 0.0)
multi_snapshot(cancer, time_steps = 9000, interval = 1000, destination = "paper_grid_dosage")

cancer = initialize(3000, 43, dosage_interval = 10, di2 = 80, dosage = 0.5, alpha = 0.0, beta = .1, gamma = 2.0, birth_population = 100, tail_skew = .5, initial_velocity = 0.1, initial_variance = 3.0, velocity = 0.1, evolution_rate = 0.0)
cancer.treatment_function = dumb_vacation_treatment!
multi_snapshot(cancer, time_steps = 9000, interval = 1000, destination = "dumb_dosage")
# cancer = initialize(3000, 42, treatment_function = decreasing_adaptive_treatment!, dosage_interval = 70)
# multi_snapshot(cancer, time_steps = 2000, interval = 300, destination = "decreasing_adaptive_treatment_multi")
# # video(cancer, time_steps = 1500, destination = "decreasing_adaptive_treatment.mp4")

# # cancer = initialize(3000, 42, treatment_function = smooth_adaptive_treatment!)
# # video(cancer, time_steps = 1500, destination = "smooth_adaptive_treatment.mp4")

# # cancer = initialize(3000, 42, dosage_interval = 1000000)
# # video(cancer, time_steps = 1500, destination = "maximum_tolerated_dosage.mp4")

cancer = initialize(3000, 44, dosage_interval = 40, birth_population = 50, mean = 1.0, std = 0.25, evolution_rate = 0.0, initial_velocity = .01)
# video(cancer, time_steps = 1500, destination = "frequent_update.mp4")
snapshot(cancer, destination = "frequent_update_snapshot.png")

# cancer = initialize(9000, 42, size = 200)
# video(cancer, break_condition = 20000, time_steps = 2000, destination = "large.mp4")
c_means = [1.0, 1.0, 1.0]
c_stds = [0.25, 0.20, 0.3]
prefixes = ["B_grid", "C_grid", "D_grid"]

evolution_rate = 0.1
velocity = 0.1

for i in 1:length(c_means)
    c_mean = c_means[i]
    c_std = c_stds[i]
    prefix = prefixes[i]
    
    model_smooth = initialize(4000, 20, alpha=0.7, beta=1.0, birth_population = 100, mean = c_mean, std = c_std, evolution_rate = evolution_rate, initial_velocity = .01, velocity = velocity, size = 150, dosage = .5, treatment_function = smooth_adaptive_treatment!)
    modelat2 = initialize(4000, 20, alpha=0.5, beta=0.1, birth_population = 100, mean = c_mean, std = c_std, evolution_rate = evolution_rate, initial_velocity = .01, velocity = velocity, size = 150, dosage = .5)
    model_slow_vacation = initialize(4000, 20, dosage_interval = 12, di2 = 120, alpha = 0.5, beta = .5, birth_population = 100, mean = c_mean, std = c_std, evolution_rate = evolution_rate, initial_velocity = .01, velocity = velocity, size = 150, dosage = .5, treatment_function = smart_vacation_treatment!)
    model_fast_vacation = initialize(4000, 20, dosage_interval = 12, di2 = 72, alpha = 0.5, beta = .5, birth_population = 100, mean = c_mean, std = c_std, evolution_rate = evolution_rate, initial_velocity = .01, velocity = velocity, size = 150, dosage = .5, treatment_function = smart_vacation_treatment!)
    modelat3 = initialize(4000, 20, alpha=0.8, beta=0.6, birth_population = 100, mean = c_mean, std = c_std, evolution_rate = evolution_rate, initial_velocity = .01, velocity = velocity, size = 150, dosage = .5)
    


    multi_snapshot(
        [model_smooth, modelat2, model_slow_vacation, model_fast_vacation, modelat3];
        destinations = ["$(prefix)_treatment_SMO", "$(prefix)_treatment_AT2", "$(prefix)_treatment_SV", "$(prefix)_treatment_FV", "$(prefix)_treatment_AT3"],
        time_steps = 19000,
        interval = 200,
        break_condition = 20000,
        combined_name = "$(prefix)_treatment_comparison"
    )
end




cancer = initialize(10000, 42, dosage_interval = 12, di2 = 72, dosage = .3, alpha = 0.5, beta = .2, gamma =  .1, birth_population = 100, mean = 1, std = 0.25, evolution_rate = .1, initial_velocity = .01, velocity = .1, size = 160, treatment_function = smart_vacation_treatment!)
cancer.dosage
cancer.last_dosage
cancer.initial_agents = 5000
video(cancer, time_steps = 1500, destination = "decreasing.mp4", break_condition = 17000, time_resolution = 10)
multi_snapshot([cancer]; time_steps = 9000, interval = 300, destinations = ["decreasing_dosage"], break_condition = 20000, combined_name = "z")
cancer.gamma