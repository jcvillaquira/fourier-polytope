using Polymake
using LinearAlgebra

include("src/fourier.jl")
include("src/polytopes.jl")

# Parameters and polytopes.
P = cube(2)
Q = polytope.polarize(P)
R = polytope.Polytope(POINTS = vertices_array([0, 0], [1, 0], [0, 1]))

# Compute Fourier transform.
FP = fourier_transform(P)
FQ = fourier_transform(Q)
FR = fourier_transform(R)
