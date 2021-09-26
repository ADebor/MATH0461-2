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

inf = zeros(M)
sup = zeros(M)
for i = 1:M
tmp = lp_rhs_perturbation_range(constraint[i])
inf[i] = tmp[1]
sup[i] = tmp[2]

println(i)
end

plot(x_axis, inf, seriestype = :scatter, title = "Allowed perturbations of the rhs coefficients")
plot_rhs = plot!(x_axis, sup, seriestype = :scatter)
savefig(plot_rhs, "L1_rhs_uncorrupted_608.png")

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

FileIO.save(File(format"PNG", "Cell_uncorrupted_608_L1.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_1014_L1.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_1521_L1.png"), img)
#FileIO.save(File(format"PNG", "Cell_uncorrupted_3042_L1.png"), img)

#= d = dual.(constraint)
val_mean_d = Statistics.mean(d)
mean_d = zeros(M)
for i = 1:M
  mean_d[i] = val_mean_d
end

plot(x_axis, d, seriestype = :scatter)
plotd = plot!(x_axis, mean_d)

savefig(plotd, "L1_dual_uncorrupted_608.png")
=#
