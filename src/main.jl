using Agents
using Random
using Distributions
using StatsBase
using ColorSchemes
using CairoMakie
using Printf
using Statistics: mean

import FromFile: @from 
@from "setup.jl" import initialize
@from "visualize.jl" import video, snapshot, multi_snapshot
@from "treatment.jl" import default_adaptive_treatment!, smooth_adaptive_treatment!, decreasing_adaptive_treatment! 

cancer = initialize(3000, 42)
# snapshot(cancer, destination = "initial_snapshot.png")
# # video(cancer, time_steps = 1500, destination = "adaptive_treatment.mp4")

cancer = initialize(3000, 42)
multi_snapshot(cancer, time_steps = 1500, interval = 300, destination = "adaptive_treatment_multi")

cancer = initialize(3000, 43, dosage_interval = 10, di2 = 80, dosage = 0.5, alpha = 0.7, beta = .2, gamma = 1.0, birth_population = 100, tail_skew = .5, initial_velocity = .1, initial_variance = 3.0, velocity = 0.1)
# cancer.treatment_function = decreasing_adaptive_treatment!
video(cancer, time_steps = 1500, destination = "decreasing_dosage_update.mp4")
multi_snapshot(cancer, time_steps = 2000, interval = 300, destination = "frequent_dosage_update_multi")

# cancer = initialize(3000, 42, treatment_function = decreasing_adaptive_treatment!, dosage_interval = 70)
# multi_snapshot(cancer, time_steps = 2000, interval = 300, destination = "decreasing_adaptive_treatment_multi")
# # video(cancer, time_steps = 1500, destination = "decreasing_adaptive_treatment.mp4")

# # cancer = initialize(3000, 42, treatment_function = smooth_adaptive_treatment!)
# # video(cancer, time_steps = 1500, destination = "smooth_adaptive_treatment.mp4")

# # cancer = initialize(3000, 42, dosage_interval = 1000000)
# # video(cancer, time_steps = 1500, destination = "maximum_tolerated_dosage.mp4")

cancer = initialize(3000, 45, dosage_interval = 40, birth_population = 200, tail_skew = 1, initial_velocity = .6, initial_variance = 10.)
# video(cancer, time_steps = 1500, destination = "frequent_update.mp4")
snapshot(cancer, destination = "frequent_update_snapshot.png")

# cancer = initialize(9000, 42, size = 200)
# video(cancer, break_condition = 20000, time_steps = 2000, destination = "large.mp4")



