include("setup.jl")
include("visualize.jl") 

cancer = initialize(3000, 42)
video(cancer, time_steps = 1500, destination = "adaptive_cancer.mp4")

cancer = initialize(3000, 42, dosage_interval = 1000000)
video(cancer, destination = "no_adaptive_cancer.mp4")

cancer = initialize(3000, 42, tail_skew = .3, dosage = 0)
video(cancer, time_steps = 50, destination = "front_cancer.mp4")