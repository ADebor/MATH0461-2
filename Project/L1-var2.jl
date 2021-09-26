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
theta = phi * psi;   # M x Nµ

lambda = 100

# CREATE EMPTY MODEL

model = Model(Gurobi.Optimizer);

# ADD VARIABLES

@variable(model, x[1:N]);
@variable(model, u[1:N]);
@variable(model, v);

# ADD CONSTRAINTS

@constraint(model, [v; theta*x - m] in SecondOrderCone())

for i = 1:N
  @constraint(model, x[i] - u[i] <= 0);
  @constraint(model, -x[i] - u[i] <= 0);
end

obj = sum(u) + lambda * v;
@objective(model, Min, obj);

# SOLVE MODEL

optimize!(model);

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

FileIO.save(File(format"PNG", "Cell_noisy_608_var2_100.png"), img)
#FileIO.save(File(format"PNG", "Cell_noisy_1014_var2_100.png"), img)
#FileIO.save(File(format"PNG", "Cell_noisy_1521_var2_100.png"), img)
#FileIO.save(File(format"PNG", "Cell_noisy_3042_var2_100.png"), img)

#ft = fftshift(fft(img))
#heatmap(log.(abs.(ft.*ft)).+1)
