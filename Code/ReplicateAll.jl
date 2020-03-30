# This is the main script that runs the simulations. 
# Written and tested in Julia 1.4.0

using Parameters
using SparseArrays
using SpecialFunctions
using LinearAlgebra
using Distributions
using QuadGK
using DelimitedFiles
using CSV
using PGFPlotsX

output_dir = joinpath("..", "Output")

include("Parameters.jl") # This contains the parameters of the model
include("Structs.jl")
include("IncomeAndGrids.jl")
include("SavingsSafeZone.jl")
include("SavingsCrisisZone.jl")
include("FindSEqms.jl")
include("FindBEqms.jl")
include("SBPolicies.jl")
include("BorrowingFull.jl")
include("BorrowingSimulation.jl")
include("PolicyAndMoments.jl")
include("MakeFigures.jl")

# Parameters
#
# The following sets the parameters of Chatterjee and Eyigungor (2012)
# plus the disaster state specification, and additional parameters used for
# the simulations (grid points, safe zone position, etc). 
par = ParNT()

# Obtaining the income process that fits the target moments
yParams, mParams, fitted_rho_eta = obtain_income_process(par)

# Looking for a savings equilibrium for the benchmark maturity (1/delta)
matListSavingsEqm, = find_savings_equilibria(
    [par.benchmark_maturity], par, yParams, mParams,
    waves=[(mixQFac=0.75, mixVFac=0.75, maxIter=4_000)] # using the smoothing parameters 
    # that we know work. 
)

# Getting the policies
if par.benchmark_maturity in matListSavingsEqm
     save_policy_moments_savings(par.benchmark_maturity, par, yParams, mParams)
     makeFigures(par, par.benchmark_maturity, "savings")
end

# Looking for borrowing equilibria for the benchmark maturity (1/delta)
matListBorrowingEqm, = find_borrowing_equilibrium(
    [par.benchmark_maturity], par, yParams, mParams,
    waves=[(mixQFac=0.98, mixVFac=0.98, maxIter=30_000)] # using the smoothing parameters 
    # that we know work. 
)

# Getting the policies
if par.benchmark_maturity in matListBorrowingEqm
    save_policy_moments_borrowing(par.benchmark_maturity, par, yParams, mParams)
    makeFigures(par, par.benchmark_maturity, "borrowing")
end

# Testing existence for other maturity levels (this can take a long time)
println("\n\n\n")
println("Looking for equilibria for many maturities. This takes a long time.")

# maturity range -- this corresponds to 1/delta
matList = [2.0, 9.0, 10.0, 30.0, 33.0, 35.0] 

#  Looking for a savings equilibrium
matListSavingsEqm, = find_savings_equilibria(
    matList, par, yParams, mParams; 
    save_output=false, 
    restart_waves=false # assuming that it gets harder for longer maturities
                        # and matList is sorted in increasing order
)

# Looking for borrowing equilibria
matListBorrowingEqm, = find_borrowing_equilibrium(
    reverse(matList), par, yParams, mParams;
    save_output=false,
    restart_waves=false  # assuming that it gets harder for shorter maturities
                         # and matList is sorted in decreasing order
)

@info "Total maturities searched" matList
@info "Multiple equilibria found for maturities in" intersect(matListSavingsEqm, matListBorrowingEqm)