# This file contains the basic parameters used for the simulations. The 
# basic parameters are obtained from Chaterjee-Eyigungor (2012) AER.

ParNT = @with_kw ( 
    beta = 0.9540232420,
    theta = 0.0385,
    gamma = 2.0,
    hpen0 = -0.1881927550,
    hpen1 = 0.2455843389,
    R = 1.01,

    benchmark_maturity = 20.0,
    delta = 1 / benchmark_maturity,  # maturity is 20 periods (quarters)
    kappa = 0.03,  # coupon 

    y_points = 201,  # value for testing
    # y_points = 201,  # value for paper's results

    # target parameters for the income process
    rho_target = 0.948503,
    eta_squared_target = 0.027092^2,

    # other parameters for the income process
    std_span = [-3.0, 3.0],
    inflate_end_points = false,
    ergodic_rentry = false,
    # m-shock specification
    m_points = 12,
    eps_squared = 0.003^2,
    m_mu = 0.0,
    m_std_span = 2.0,
    # getting the transition matrix for the disaster state
    pi_LL = 1.0 - 1.0 / (3.5 * 4.0),
    pi_HH = (1.0 - 0.0383)^(1/4),
    mu_T = [pi_LL (1.0 - pi_HH); (1.0 - pi_LL) pi_HH],
    
    # A disaster drops output by an average of 0.20
    # and the unconditional mean of log output is 0.0
    mu_system_soln = [0.20, 0.0],  

    # Additional parameters
    aLim = -0.89,  # Safe Zone bound
    aExtBoundMag = 1.0, # extend the crisis zone to -1.89
    a_crisis_zone_points = 500,
    a_total_points = 946,

    # Parameters for the grids derived from above
    aPoints = a_total_points - a_crisis_zone_points,  # points in the Safe Zone
    aExtPoints = a_crisis_zone_points,  # points in the Crisis Zone
    aBoundsSafe = [aLim, 0.0],  # Safe Zone region
    aFinePct = 0.9,  # There will be a finer grid 10% to the right of aLim .. 
    aFinePts = 89,   # .. with this number of points. 

    # directories to save output files 
    resultDirSave = joinpath(output_dir, "Models", "Savings"),
    resultDirBorrow = joinpath(output_dir, "Models", "Borrowing"),
    resultDirFigures = joinpath(output_dir, "Figures"),
    resultDirMoments = joinpath(output_dir, "Moments"),
    resultDirPolicies = joinpath(output_dir, "Policies")
)


# u is the CRRA utility function. Its arguments are just:
# 1. m: the model specification
# 2. x: the value of consumption
function u(m, x)
    # If x is positive return the value of CRRA utility
    if x > 0.0
        if m.gamma != one(m.gamma)
            return x^(1.0 - m.gamma) / (1.0 - m.gamma)
        else
            return log(x)
        end
    # Otherwise return either a very negative value, scaled by how negative x is when utility is negative near x=0
    # or a large negative value which gets smaller as x rises to 0.
    # This is done to ensure that any optimization algorithm knows to try to make consumption positive at essentially any cost.
    else
        if m.gamma >= one(m.gamma)
            return u(m, 1e-10) * (1.0 + abs(x))
        else
            return -1e10 * abs(x)
        end
    end
end

# uInverse is the inverse of the CRRA utility function. Its arguments are just:
# 1. m: the model specification
# 2. x: the value of consumption

# UNLIKE THE FUNCTION JUST ABOVE, THIS FUNCTION ASSUMES THAT gamma > 1.
function uInverse(m, x)
    if x < 0.0
        return ((1.0 - m.gamma) * x)^(1.0 / (1.0 - m.gamma))
    else
        return uInverse(m, -1e-10) / (1.0 + abs(x))
    end
end

