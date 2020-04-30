# Self-fulfilling Debt Dilution Repository

This repository contains the source code for replicating the results in the paper:

   Aguiar, Mark and Manuel Amador (2020): “Self-fulfilling Debt Dilution: 
   Maturity and Multiplicity in Debt Models”.

The authors (Mark Aguiar and Manuel Amador) are grateful for the fantastic research assistance of 
Stelios Fourakis in writing and testing this code. All remaining errors are the responsibility 
of the authors.

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

The code is in Julia: https://julialang.org/

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

https://docs.julialang.org/en/v1/manual/parallel-computing/


_WARNING:_ For the figures, you need a working latex installation (see 
https://kristofferc.github.io/PGFPlotsX.jl/stable/).  In case you do not 
have one, comment out the calls to `makeFigures` in the `ReplicateAll.jl` 
script in lines 57 and 70 before running the script.

### To run the script

Open a terminal. Make sure that your working directory is `Code`. If not, change
the working directory to the `Code` subdirectory. Start julia and use the 
following command at the julia REPL prompt:
    
    julia> include("ReplicateAll.jl")

## Other relevant files 

`Parameters.jl` contains the parameters used for the simulations. 

The figures shown in the paper are saved in the `Output/Figures` directory. 

The moments used to construct the tables in the paper are saved in the 
`Output/Moments` directory. Each moments file contains two columns. 
The first represents the moments using the full sample. 
The second column represents the moments using excluding the disaster states. 
The moments are:

    1. mean A' over Y no default: E[A'/(y+m)] when the country has been out of default long 
       enough and does not default today.

    2. mean market value over Y: E[q(y,A')*A'/(y+m)] when the country has been out of default 
       long enough and does not default today.

    3. mean A over Y no default: E[A/(y+m)] when the country has been out of default long 
       enough and does not default today.

    4. mean A over Y: E[A/(y+m)] when the country has been out of default 
       long enough up until today.

    5. def Rate: E[d] when the country has been out of default long enough 
       up until today.

    6. mean Spread: E[(1+r(y,A'))^4-(1+r)^4] when the country has been out of 
       default long enough and does not default today.

    7. vol Spread: SD[(1+r(y,A'))^4-(1+r)^4] when the country has been out of 
       default long enough and does not default today.

    8. vol C over Y: SD[ln(c)]/SD[ln(y+m)] when the country has been out of 
       default long enough and does not default today.

    9. vol TB: SD[(y+m-c)/(y+m)] when the country has been out of default long 
       enough and does not default today.

    10. cor TB with Log Y: corr[(y+m-c)/(y+m),ln(y+m)] when the country has been out of 
        default long enough and does not default today.

    11. cor Spread with Log Y: corr[(1+r(y,A'))^4-(1+r)^4,ln(y+m)] when the country has 
        been out of default long enough and does not default today.

    12. cor Spread with A / Y: corr[(1+r(y,A'))^4-(1+r)^4,A'/(y+m)] when the country 
        has been out of default long enough and does not default today.

    13. cor Spread with TB: corr[(1+r(y,A'))^4-(1+r)^4, (y+m-c)/(y+m) ] when the 
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
additional computations with 6 more maturities run just under 17h.

