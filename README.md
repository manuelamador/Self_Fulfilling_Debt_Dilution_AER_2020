# Self-fulfilling Debt Dilution Repository

This repository contains the source code for replicating the results in the paper:

   Aguiar, Mark and Manuel Amador (2020): “Self-fulfilling Debt Dilution: 
   Maturity and Multiplicity in Debt Models”.

The authors (Mark Aguiar and Manuel Amador) are grateful for the fantastic research assistance of 
Stelios Fourakis in writing and testing this code. All remaining errors are the responsibility 
of the authors.

## Summary 

This repository contains the code to generate the numerical simulation results
presented in Section 8 and Appendix E. The present code generates the four plots shown 
in the two figures contained in Appendix E as well as the corresponding table of moments. 

_NOTE:_ The other figures in the paper are not the result of this numerical
simulation, and were directly obtained from closed-form formulas using the
parameter values specified in the caption for each figure. 

An up-to-date version of this repository can be found in: 

[https://github.com/manuelamador/Self_Fulfilling_Debt_Dilution_AER_2020](https://github.com/manuelamador/Self_Fulfilling_Debt_Dilution_AER_2020)


## Folder structure 

The folder structure is as follows:

- `Code`: contains all of the julia source code necessary for replication
- `Output`: contains the output files from the simulations
    - `Figures`: contains the generated figures
    - `Moments`: contains the moments of each simulation
    - `Policies`: contains the policies for each simulation
    - `Models/Savings`: contains the model elements for a savings simulation
    - `Models/Borrowing`: contains the model elements for a borrowing simulation

## Running the code

The code is in Julia: [https://julialang.org/](https://julialang.org/)

The main source file is `ReplicateAll.jl`. This script runs the simulations for
both the borrowing and the saving equilibrium and generates the moments
and figures. 

_WARNING:_ You'll need the following Julia packages installed: 

    Parameters,
    SparseArrays,
    SpecialFunctions,
    LinearAlgebra,
    Distributions,
    QuadGK,
    DelimitedFiles,
    CSV,
    PGFPlotsX

_WARNING:_ Before running this script, make sure you have the appropriate 
directory structure shown above. 

_WARNING:_ Make sure that you allow Julia to use multiple threads: 

[https://docs.julialang.org/en/v1/manual/parallel-computing/](https://docs.julialang.org/en/v1/manual/parallel-computing/)


_WARNING:_ For the figures, you need a working latex installation (see 
[https://kristofferc.github.io/PGFPlotsX.jl/stable/](https://kristofferc.github.io/PGFPlotsX.jl/stable/)). 
In case you do not 
have one, comment out the calls to `makeFigures` in the `ReplicateAll.jl` 
script in lines 57 and 70 before running the script.

### To run the script

Open a terminal. Make sure that your working directory is `Code`. If not, change
the working directory to the `Code` subdirectory. Start julia and use the 
following command at the julia REPL prompt:
    
    julia> include("ReplicateAll.jl")

## Other relevant files 

### Parameters

`Parameters.jl` contains the parameters used for the simulations. 

### Figures 

The figures shown in the paper are saved in the `Output/Figures` directory:

   1. `pol_S_plot.pdf`: Corresponds to panel (a) of the figure "Simulation Results: Policy Functions"
   2. `pol_B_plot.pdf`: Corresponds to panel (b) of the figure "Simulation Results: Policy Functions"
   3. `price_S_plot.pdf` : Corresponds to panel (a) of the figure "Simulation Results: Price Functions"
   4. `price_B_plot.pdf`: Corresponds to panel (b) of the figure "Simulation Results: Price Functions"

The function `makeFigures` defined in `makeFigures.jl` generates these figures using as an 
input the computation results stored in `Output/Models` and `Output/Policies`.

### Moments

The moments used to construct the table in the paper are saved in the 
`Output/Moments` directory. The corresponding moment files are:

   1. `mat_20.0_savings_moments.csv`: Contains the moments for the savings equilibrium.
   2. `mat_20.0_borrowing_moments.csv`: Contains the moments for the borrowing equilibrium. 

The functions `save_policy_moments_savings` and `save_policy_moments_borrowing`
defined in `PolicyAndMoments.jl` use previous computations to simulate the
model and compute and save the moments. 

Each moment file contains two columns. The first represents the moments using
the full sample. The second column represents the moments using excluding the
disaster states. The moments are:

  1. _mean A' over Y no default_: `E[A'/(y+m)]` when the country has been out of default long 
     enough and does not default today.

  2. _mean market value over Y_: `E[q(y,A')*A'/(y+m)]` when the country has been out of default 
     long enough and does not default today.

  3. _mean A over Y no default_: `E[A/(y+m)]` when the country has been out of default long 
     enough and does not default today.

  4. _mean A over Y_: `E[A/(y+m)]` when the country has been out of default 
     long enough up until today.

  5. _def Rate_: `E[d]` when the country has been out of default long enough 
     up until today.

  6. _mean Spread_: `E[(1+r(y,A'))^4-(1+r)^4]` when the country has been out of 
     default long enough and does not default today.

  7. _vol Spread_: `SD[(1+r(y,A'))^4-(1+r)^4]` when the country has been out of 
     default long enough and does not default today.

  8. _vol C over Y_: `SD[ln(c)]/V[ln(y+m)]` when the country has been out of 
     default long enough and does not default today.

  9. _vol TB_: `SD[(y+m-c)/(y+m)]` when the country has been out of default long 
     enough and does not default today.

  10. _cor TB with Log Y_: `corr[(y+m-c)/(y+m),ln(y+m)]` when the country has been out of 
      default long enough and does not default today.

  11. _cor Spread with Log Y_: `corr[(1+r(y,A'))^4-(1+r)^4,ln(y+m)]` when the country has 
      been out of default long enough and does not default today.

  12. _cor Spread with A / Y_: `corr[(1+r(y,A'))^4-(1+r)^4,A'/(y+m)]` when the country 
      has been out of default long enough and does not default today.

  13. _cor Spread with TB_: `corr[(1+r(y,A'))^4-(1+r)^4, (y+m-c)/(y+m)]` when the 
      country has been out of default long enough and does not default today.

NOTE: the code is written where `A` represents assets. The paper instead uses 
the convention and writes the relevant moments and figures using debt, `b = -a`. 
So the relevant moments need to change sign. Note also that `A'` represents
end of period assets.


## Timing

The code was tested in Amazon Web Services (AWS) using a 
c4.4xlarge linux instance with 16 virtual CPUs and approximately 32GB of RAM.

The main part of the code (the computation of the equilibria for the benchmark
maturity of 1/20.0) ran under 1h40minutes. The entire code including the 
additional computations with additional maturities ran just under 17h.

