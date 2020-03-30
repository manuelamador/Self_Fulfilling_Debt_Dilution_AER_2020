function find_borrowing_equilibrium(
    invMatParamList,
    par,
    yParams, 
    mParams;  
    save_output=true,
    restart_waves::Bool=true,
    verbose::Bool=true,
    waves=[
        (mixQFac=0.95, mixVFac=0.95, maxIter=15_000),
        (mixQFac=0.98, mixVFac=0.98, maxIter=30_000),
        (mixQFac=0.99, mixVFac=0.99, maxIter=80_000)
    ]
)
    @unpack (
        aLim,
        aPoints,
        aFinePts,
        aFinePct,
        aExtPoints,
        aExtBoundMag,
        aBoundsSafe
    ) = par

    resultDir = par.resultDirBorrow

    mBase = LongTermBondRDSaveSpec(
        par.beta,
        par.gamma,
        par.R,
        par.delta,
        par.kappa,
        aBoundsSafe,
        aPoints,
        yParams,
        mParams
    )
    if verbose
        println("Looking for borrowing equilibria")
    end 

    matDim = length(invMatParamList)
    lambdaList = invMatParamList.^ - 1

    # risk free price
    origQBase = (
        (mBase.lambda + (1.0 - mBase.lambda) * mBase.coup ) / 
        (mBase.R - (1.0 - mBase.lambda))
    )
    modCoupList = (
        ((mBase.R .- (1.0 .- lambdaList)) * origQBase .- lambdaList) ./ 
        (1.0 .- lambdaList)
    )

    isAnEqm = zeros(Bool, matDim)
    validBEqm = zeros(Bool, matDim)
  
    aFineUB = aFinePct * aLim
    aBounds = [aLim, mBase.aBounds[2]]
    
    mSave = LongTermBondRDSaveSpec(
        mBase.beta,
        mBase.gamma,
        mBase.R,
        mBase.lambda,
        mBase.coup,
        aBounds,
        aPoints,
        mBase.yParams,
        mBase.mParams
    )
    s = makeEval(mSave, aFineUB, aFinePts)

    vfiGOneStep!(mSave, s, 4e-15, 4000, 0.5)

    mat_index = 1
    total_waves = length(waves)
    current_wave_ind = 1
    
    verbose && println()

    while true 
        maturity = invMatParamList[mat_index]
        tempLambda = lambdaList[mat_index]
        tempCoup = modCoupList[mat_index]

        next_mat = false 

        if verbose
            println("    Looking for borrowing eqm for mat ", 
                maturity, " using smoothing ", waves[current_wave_ind])
        end
        
        mMod = LongTermBondRDSaveSpec(
            mBase.beta,
            mBase.gamma,
            mBase.R,
            tempLambda,
            tempCoup,
            aBounds,
            aPoints,
            mBase.yParams,
            mBase.mParams
        )

        tempAExtBounds = [aBounds[1] - aExtBoundMag, aBounds[1]]

        sExt = setupEqmCompletionYDefRDImpl(
            mMod,
            s,
            tempAExtBounds,
            aExtPoints,
            par.theta,
            par.hpen0,
            par.hpen1
        )

        aBoundsB = [tempAExtBounds[1], aBounds[2]]
        aPointsB = aExtPoints + aPoints

        mModB = LongTermBondRDBorrowSpec(
            mBase.beta,
            par.theta,
            mBase.gamma,
            mBase.R,
            tempLambda,
            tempCoup,
            par.hpen0,
            par.hpen1,
            aBoundsB,
            aPointsB,
            mBase.yParams,
            mBase.mParams,
            false,
            waves[current_wave_ind].mixQFac,
            waves[current_wave_ind].mixVFac
        )
        sBorrow = makeEval(
            mModB,
            mMod,
            s,
            sExt
        )
        sBorrow, g2, tempMaxDist = vfiGOneStep!(
            mModB, 
            sBorrow, 
            1e-8, 
            waves[current_wave_ind].maxIter
        )
        # computing the distance with respect to the risk free price
        diffWithRFAt0 = abs.(sBorrow.qGrid[sBorrow.a0Ind, :] .- origQBase)
        verbose && println()
        
        @debug "parameters" maturity aLim tempMaxDist maximum(diffWithRFAt0)

        if (tempMaxDist < (1e-6))
            # found an equilibrium
            verbose && println("    Found an equilibrium... ")
            if (maximum(diffWithRFAt0) >= (1e-4))
                validBEqm[mat_index] = true            
                next_mat = true
                verbose && println("    and it is a borrowing one. Saving to disk.") 
            else
                # But it is not a valid borrowing one (it is a savings one!) as
                # the price at zero is the risk-free price. 
                verbose && println("    but IT IS NOT a borrowing one.")
            end
            paramListOut = [
                invMatParamList[mat_index],
                lambdaList[mat_index],
                modCoupList[mat_index],
                aLim,
                tempMaxDist,
                maximum(diffWithRFAt0),
                validBEqm[mat_index]
            ]
            if save_output && validBEqm[mat_index]
                writecsv(
                    joinpath(resultDir, "mat_$maturity"*"_qAt0BAndLim.csv"),
                    hcat(
                        sBorrow.qGrid[aExtPoints + 1, :],
                        sBorrow.qGrid[sBorrow.a0Ind, :]
                    )
                )
                writecsv(
                    joinpath(resultDir, "mat_$maturity"*"_params.csv"),
                    paramListOut
                )
                writeEssentials(sBorrow, "mat_$maturity", resultDir)
            end
        end 
        if !next_mat 
            if current_wave_ind < total_waves
                current_wave_ind += 1 
                # increase current wave 
                verbose && println("    Didn't find eqm. Trying higher smoothing pars.")
            else
                verbose && println("    Didn't find eqm with highest smoothing pars.")
                next_mat = true
            end 
        end   
        if next_mat
            (mat_index >= matDim) && break
            # moving on the next maturity index
            mat_index = mat_index + 1   
            # restarting the wave of smoothing parameters if required
            restart_waves && (current_wave_inds = 1)
        end
    end

    noBorrowinEqmMaturities = invMatParamList[.!validBEqm]
    matListBorrowingEqm = invMatParamList[validBEqm]

    if verbose
        @info "Done. Found borrowing equilibria for $matListBorrowingEqm"
        if !isempty(noBorrowinEqmMaturities) 
            @warn "But couldn't compute borrowing eqm for $noBorrowinEqmMaturities"
        end
    end

    return matListBorrowingEqm, noBorrowinEqmMaturities
end

