using JuMP
using Gurobi
using MosekTools
using MathOptInterface
using ImageView
using Images
using FileIO
using FFTW
using Plots

include("utilities.jl")

# DEFINE PARAMETERS

m = unpickler("uncorrupted_measurements_M608.pickle")        # M
#m = unpickler("uncorrupted_measurements_M1014.pickle")
#m = unpickler("uncorrupted_measurements_M1521.pickle")
#m = unpickler("uncorrupted_measurements_M3042.pickle")
phi = unpickler("measurement_matrix_M608.pickle")           # M x N
#phi = unpickler("measurement_matrix_M1014.pickle")
#phi = unpickler("measurement_matrix_M1521.pickle")
#phi = unpickler("measurement_matrix_M3042.pickle")
psi = unpickler("basis_matrix.pickle")                      # N x N

M = size(m,1)
N = size(phi, 2)
theta = phi * psi   # M x N

# CREATE EMPTY MODEL

LP_model = Model(Gurobi.Optimizer)

# ADD VARIABLES

@variable(LP_model, x[1:N])

# ADD CONSTRAINTS

@constraint(LP_model, theta*x - m .== 0)

# ADD OBJECTIVE

obj = count(i->(i!=0),x)
@objective(LP_model, Min, obj)

# WRITE MODEL TO FILE

#write_to_file(LP_model, "L0.mps")

# SOLVE MODEL

optimize!(LP_model)

# SHOW IMAGE AND SPECTRUM

x = value.(x)
r_prime = psi * x;

minr = minimum(r_prime)
for i = 1:N
    r_prime[i] = (r_prime[i]-minr)
end
maxr = maximum(r_prime)
for i = 1:N
    r_prime[i] = r_prime[i]/maxr
end

img = reshape(r_prime,(78,78));

imshow(img)

FileIO.save(File(format"PNG", "Cell_uncorrupted_608_L0.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_1014_L0.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_1521_L0.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_3042_L0.png"), img)
