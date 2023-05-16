using Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets

# We don't use BenchmarkTools for this model

model = "Beam"
cases = ["CB22d-100-C", # constant input set
         "CB22d-100-F", # bounded but arbitrarily-varying input set
         "CB22d-100-C-discrete",
         "CB22d-100-F-discrete",
         "CB22d-500-C",
#          "CB22d-500-F",
         "CB22d-500-C-discrete",
         "CB22d-500-F-discrete",
         "CB22d-1000-C",
#          "CB22d-1000-F",
#          "CB22d-1000-C-discrete",
#          "CB22d-1000-F-discrete"
        ]
if TEST_LONG
    insert!(cases, 6, "CB22d-500-F")
    push!(cases, "CB22d-1000-C-discrete")
    push!(cases, "CB22d-1000-F-discrete")
end

# measure = [] # maximum value of the velocity at node 70
validation = ones(Int, length(cases)) # nothing to verify
results = Dict(model => Dict(c => -1.0 for c in cases))
sol = Dict{String,Any}()

include("clamped.jl")
include("clamped_krylov.jl")

LazySets.deactivate_assertions()

# ----------------------------------------
#  CB22d-100-C
# ----------------------------------------

case = "CB22d-100-C"
m, N, Z = parse_clamped_beam(case)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-6, vars=[70, 170], n=201)
res = @timed solve_krylov_discr(prob, NSTEPS=10_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CB22d-100-C-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[70, 170], n=201, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CB22d-100-F
# ----------------------------------------

case = "CB22d-100-F"
m, N, Z = parse_clamped_beam(case)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-7, vars=[70, 170], n=200)
res = @timed solve(prob, NSTEPS=100_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CB22d-100-F-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[70, 170], n=200, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CB22d-500-C
# ----------------------------------------

case = "CB22d-500-C"
m, N, Z = parse_clamped_beam(case)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-6, vars=[350, 850], n=1001)
res = @timed solve_krylov_discr(prob, NSTEPS=10_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CB22d-500-C-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[350, 850], n=1001, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CB22d-500-F
# ----------------------------------------

case = "CB22d-500-F"
m, N, Z = parse_clamped_beam(case)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-7, vars=[350, 850], n=1000)
res = @timed solve(prob, NSTEPS=100_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CB22d-500-F-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[350, 850], n=1000, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CB22d-1000-C
# ----------------------------------------

case = "CB22d-1000-C"
m, N, Z = parse_clamped_beam(case)
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-6, vars=[700, 1700], n=2001)
res = @timed solve_krylov_discr(prob, NSTEPS=10_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

case = "CB22d-1000-C-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[700, 1700], n=2001, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ----------------------------------------
#  CB22d-1000-F
# ----------------------------------------

case = "CB22d-1000-F"
m, N, Z = parse_clamped_beam(case)
#=
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=1e-7, vars=[700, 1700], n=2000)
res = @timed solve(prob, NSTEPS=100_000, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time
=#

case = "CB22d-1000-F-discrete"
prob = clamped(case, mode=m, N=N, Z=Z)
alg = LGG09(δ=9.88e-7, vars=[700, 1700], n=2000, approx_model=NoBloating())
res = @timed solve(prob, T=0.01, alg=alg)
sol[case] = res.value
# push!(measure, ρ(alg.template.directions[2], sol[case]))
results[model][case] = res.time

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# export runtimes
runtimes = Dict()
for (i, c) in enumerate(cases)
   local t = results[model][c] # in seconds
   runtimes[c] = round(t, digits=4)
end

if !@isdefined io
    io = stdout
end

for (i, c) in enumerate(cases)
   print(io, "$model,$c,$(validation[i]),$(runtimes[c])\n")
end

# ==============================================================================
# Plot
# ==============================================================================

if !@isdefined TARGET_FOLDER
    TARGET_FOLDER = @__DIR__
end

# Constant force case
fig = Plots.plot()
Plots.plot!(fig, sol["CB22d-100-C"], vars=(0, 170),
           color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
           tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"t",
           ylab=L"v_{70}",
           xtick=([0, 0.0025, 0.005, 0.0075, 0.01],
                  [L"0.0", L"0.0025", L"0.0050", L"0.0075", L"0.0100"]),
           ytick=([-80, -60, -40, -20, 0, 20, 40, 60, 80],
                  [L"-80", L"-60", L"-40", L"-20", L"0", L"20", L"40", L"60", L"80"]),
           xlims=(0, 0.01), ylims=(-80, 80),
           bottom_margin=0mm, left_margin=0mm, right_margin=15mm, top_margin=3mm,
           size=(1000, 1000))

# savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP23-JuliaReach-$model-CB22d100C.pdf"))
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP23-JuliaReach-$model-CB22d100C.png"))

# Time-varying force case
fig = Plots.plot()
Plots.plot!(fig, sol["CB22d-100-F"], vars=(0, 170),
           color=:blue, alpha=0.5, lw=1.0, linecolor=:blue,
           tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"t",
           ylab=L"v_{70}",
           xtick=([0, 0.0025, 0.005, 0.0075, 0.01],
                  [L"0.0", L"0.0025", L"0.0050", L"0.0075", L"0.0100"]),
           ytick=([-80, -60, -40, -20, 0, 20, 40, 60, 80],
                  [L"-80", L"-60", L"-40", L"-20", L"0", L"20", L"40", L"60", L"80"]),
           xlims=(0, 0.01), ylims=(-80, 80),
           bottom_margin=0mm, left_margin=0mm, right_margin=15mm, top_margin=3mm,
           size=(1000, 1000))

# savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP23-JuliaReach-$model-CB22d100F.pdf"))
savefig(fig, joinpath(TARGET_FOLDER, "ARCH-COMP23-JuliaReach-$model-CB22d100F.png"))

sol = nothing
GC.gc()
