include("setup.jl")
include("visualize.jl") 

cancer = initialize(2000, 42)
video(cancer)

cancer = initialize(2000, 43, dosage_interval = 1000000)
video(cancer, destination = "no_adaptive_cancer.mp4")