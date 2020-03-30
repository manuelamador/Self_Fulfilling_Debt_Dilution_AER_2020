function makeFigures(par, maturity::Real, eqm_type::String; verbose=true)

    @assert eqm_type in ["savings", "borrowing"]

    if verbose
        println("Generating the figures for " * eqm_type * " eqm.")
    end

    if eqm_type == "savings"
        baselineEqm = par.resultDirSave
        underscore = "S"
    else
        baselineEqm = par.resultDirBorrow
        underscore = "B"
    end 

    baselinePolicies = par.resultDirPolicies

    aGrid=CSV.read(
        joinpath(baselineEqm, "mat_$(maturity)_aGrid.csv"), 
        header=false
    )
    yGrid=CSV.read(
        joinpath(baselineEqm, "mat_$(maturity)_yGrid.csv"), 
        header=false
    )
    q=CSV.read(
        joinpath(baselineEqm, "mat_$(maturity)_qGrid.csv"), 
        header=false
    )
    policy=CSV.read(
        joinpath(baselinePolicies, "mat_$(maturity)_$(eqm_type)_apRGrid.csv"), 
        header=false
    )
    default=CSV.read(
        joinpath(baselinePolicies, "mat_$(maturity)_$(eqm_type)_defaultProb.csv"),
        header=false
    )

    x = -aGrid[!, 1]
    risk_free = 1.0 # 1.308 

    function pol_and_def(pol, def)
        out = zero(pol)
        for i in eachindex(def)
            if def[i] < 0.1
                out[i] = pol[i]
            end
        end
        return out
    end

    max_default_risk = mapslices(
        maximum, convert(Matrix, default), dims=2)[:, 1]
    safezone = findfirst(x -> x < 0.000001, max_default_risk)

    midpoint_state = (size(policy)[2] - 1) ÷ 2 # getting the middle of the income
    disaster_state = size(policy)[2] # the disaster state is the last one
    
    lambda = 1 / maturity
    
    risk_free = (
        (lambda + (1.0 - lambda) * par.kappa ) / 
        (par.R - (1.0 - lambda))
    )

    # prices
    p1_coord1 = Coordinates(x, q[!, midpoint_state])
    p1_coord2 = Coordinates(x, q[!, end])
    p1_coord3 = Coordinates(x, risk_free * ones(length(x)))
    p1 = @pgf TikzPicture( 
        Axis(
            {
                xlabel="b'", ylabel="\\textrm{Price}", 
                legend_pos="south west" , xmin=0.0, xmax=1.4, ymin=0.0
            },
            PlotInc(
                {
                    style="black, ultra thick", mark="none"
                },
                p1_coord1
            ),
            LegendEntry("\$q_$(underscore)(b', \\mu)\$"), 
            PlotInc(
                {
                    style="black, dashed, dash pattern=on 5pt off 5pt, ultra thick", 
                    mark="none"  
                },
                p1_coord2,
            ),
            LegendEntry("\$q_$(underscore)(b', y_{dis})\$"),
            PlotInc(
                {
                    style="gray", mark="none"    
                },
                p1_coord3, 
            ),
            [
                "\\fill[blue, opacity=.2] (-1,-1) rectangle ($(x[safezone]), 2);"
            ],
            [
                "\\node ",
                " at ",
                Coordinate(1.308, 1.238),
                "{\$q^\\star\$};"
            ]
        )
    )

    # policies
    p2_coord1 = Coordinates(
        x, 
        pol_and_def(
            -policy[!, midpoint_state], 
            default[!, midpoint_state]
        )
    )
    p2_coord2 = Coordinates(
        x, 
        pol_and_def(
            -policy[!, disaster_state], 
            default[!, disaster_state]
        )
    )
    p2 = @pgf TikzPicture( 
        Axis(
            {
                xlabel="b", ylabel="b'", 
                legend_pos="north west" , xmin=0.0, xmax=1.4, ymin=0.0,
                ymax=1.2
            },
            PlotInc(
                {
                    style="black, ultra thick", mark="none"
                },
                p2_coord1
            ),
            LegendEntry("\$B_$(underscore)(b, \\mu)\$"), 
            PlotInc(
                {
                    style="black, dashed, dash pattern=on 5pt off 5pt, ultra thick", 
                    mark="none"  
                },
                p2_coord2,
            ),
            LegendEntry("\$B_$(underscore)(b, y_{dis})\$"),
            PlotInc(
                {
                    style="gray",
                    mark="none"
                },
                Coordinates(x, x)
            ),
            [
                "\\fill [blue, opacity=.2] (-1,-1) rectangle ($(x[safezone]), 2);",
                "\\node at ", 
                Coordinate(0.2, 0.13),
                "{\$45^\\circ\$};"
            ]
        )
    )

    pgfsave(joinpath(par.resultDirFigures, "price_$(underscore)_plot.pdf"), p1)
    pgfsave(joinpath(par.resultDirFigures, "pol_$(underscore)_plot.pdf"), p2)

    if verbose
        println("   Figures saved.")
    end    
end

