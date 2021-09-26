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
theta = phi * psi;   # M x N

# CREATE EMPTY MODEL

SOC_model = Model(Gurobi.Optimizer);

# ADD VARIABLES

@variable(SOC_model, x[1:N]);
@variable(SOC_model, t);

# ADD CONSTRAINTS

@constraint(SOC_model, theta*x - m .== 0);
@constraint(SOC_model,[t; x] in SecondOrderCone())

# ADD OBJECTIVE

obj = t
@objective(SOC_model, Min, obj);

# SOLVE MODEL

optimize!(SOC_model);

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

FileIO.save(File(format"PNG", "Cell_uncorrupted_608_L2.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_1014_L2.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_1521_L2.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_3042_L2.png"), img)
