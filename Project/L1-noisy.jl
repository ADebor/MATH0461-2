using JuMP
using Gurobi
using MosekTools
using MathOptInterface
using ImageView
using Images
using FileIO
using FFTW
using Plots
using Statistics

include("utilities.jl")

# DEFINE PARAMETERS

m = unpickler("noisy_measurements_M608.pickle")        # M
#m = unpickler("noisy_measurements_M1014.pickle")
#m = unpickler("noisy_measurements_M1521.pickle")
#m = unpickler("noisy_measurements_M3042.pickle")
phi = unpickler("measurement_matrix_M608.pickle")      # M x N
#phi = unpickler("measurement_matrix_M1014.pickle")
#phi = unpickler("measurement_matrix_M1521.pickle")
#phi = unpickler("measurement_matrix_M3042.pickle")
psi = unpickler("basis_matrix.pickle")                 # N x N

M = size(m,1)
N = size(phi, 2)
theta = phi * psi;   # M x N
x_axis = 1:M

# CREATE EMPTY MODEL

LP_model = Model(Gurobi.Optimizer);

# ADD VARIABLES

@variable(LP_model, x[1:N]);
@variable(LP_model, t[1:N]);

# ADD CONSTRAINTS

@constraint(LP_model, constraint, theta*x - m .== 0);

for i = 1:N
  @constraint(LP_model, x[i] - t[i] <= 0);
  @constraint(LP_model, -x[i] - t[i] <= 0);
end

# ADD OBJECTIVE

obj = sum(t)
@objective(LP_model, Min, obj);

# SOLVE MODEL

optimize!(LP_model);

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

#imshow(img)

FileIO.save(File(format"PNG", "Cell_noisy_608_L1.png"), img)
#FileIO.save(File(format"PNG", "Cell_noisy_1014_L1.png"), img)
#FileIO.save(File(format"PNG", "Cell_noisy_1521_L1.png"), img)
#FileIO.save(File(format"PNG", "Cell_noisy_3042_L1.png"), img)
