using Agents
using Random
using Distributions
using StatsBase
using ColorSchemes
using CairoMakie
using Printf

import FromFile: @from 
@from "setup.jl" import initialize
@from "visualize.jl" import video 

# cancer = initialize(3000, 42)
# video(cancer, time_steps = 1500, destination = "adaptive_treatment.mp4")

# cancer = initialize(3000, 42, treatment_function = smooth_adaptive_treatment!)
# video(cancer, time_steps = 1500, destination = "smooth_adaptive_treatment.mp4")

# cancer = initialize(3000, 42, dosage_interval = 1000000)
# video(cancer, time_steps = 1500, destination = "maximum_tolerated_dosage.mp4")

cancer = initialize(3000, 42, dosage_interval = 40)
video(cancer, time_steps = 200, destination = "frequent_update.mp4")

# cancer = initialize(9000, 42, size = 200)
# video(cancer, break_condition = 20000, time_steps = 2000, destination = "large.mp4")

cancer = initialize(3000, 42)
step!(cancer, 1000)



