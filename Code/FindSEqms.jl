function find_saving_eqm(
    mBase::LongTermBondRDSaveSpec{F,S},
    invMatParamList::Array{F,1},
    aFineUBPct::F,
    aLimit::F,
    aLimitPoint::S,
    aFinePoint::S,
    aExtPoints::S,
    aExtMag::F,
    theta::F,
    hpen0::F,
    hpen1::F,
    resultDir::String;
    waves,
    save_output::Bool=true,
    restart_waves::Bool=true,
    verbose::Bool=true
) where{F<:Real,S<:Integer} 

    lastSaveEqmExistMatInd = 1
    numNoEqmExistences = 0
    matDim = length(invMatParamList)
    lambdaList = invMatParamList.^-1
    origQBase = (mBase.lambda+(1.0 - mBase.lambda) * mBase.coup) / 
        (mBase.R - (1.0 - mBase.lambda))
    modCoupList = ((mBase.R .- (1.0 .- lambdaList)) * origQBase .- lambdaList) ./
        (1.0 .- lambdaList)

    validSaveEqm = zeros(Bool, matDim)

    alimInd  =  1
    tempALim = aLimit
    tempAFineUB = aFineUBPct*tempALim
    tempABounds = [tempALim,mBase.aBounds[2]]
    tempAPoints = aLimitPoint
    tempAFinePoints = aFinePoint
    mSave = LongTermBondRDSaveSpec(
        mBase.beta,
        mBase.gamma,
        mBase.R,
        mBase.lambda,
        mBase.coup,
        tempABounds,
        tempAPoints,
        mBase.yParams,
        mBase.mParams
    )
    s = makeEval(mSave, tempAFineUB, tempAFinePoints)
    vfiGOneStep!(mSave, s, 4e-15, 4000, 0.5)
 
    matInd = 1
    total_waves = length(waves)
    current_wave_ind = 1
    next_mat = false

    verbose && println()
    while true
        next_mat = false
        maturity = invMatParamList[matInd]
        if verbose
            println("    Looking for saving eqm for mat ", 
                maturity, " using smoothing ", waves[current_wave_ind])
        end
        tempLambda = lambdaList[matInd]
        tempCoup = modCoupList[matInd]
        current_wave = waves[current_wave_ind]

        mMod = LongTermBondRDSaveSpec(
            mBase.beta,
            mBase.gamma,
            mBase.R,
            tempLambda,
            tempCoup,
            tempABounds,
            tempAPoints,
            mBase.yParams,
            mBase.mParams
        )
        tempAExtBounds = [tempABounds[1] - aExtMag, tempABounds[1]]
        sExt = setupEqmCompletionYDefRDImpl(
            mMod,
            s,
            tempAExtBounds,
            aExtPoints,
            theta,
            hpen0,
            hpen1
        )
        sExt, g2, tempNM, tempMaxDist = vfiGOneStep!(
            mMod, s, sExt, 1e-8, 
            current_wave.maxIter, 
            current_wave.mixQFac, 
            current_wave.mixVFac
        )
        # tempNM stores whether the value function is not monotonic 
        # at the boundary of the safe zone (hence we don't have an equilibrium)
        if (tempNM == true) || (tempMaxDist >= (1e-6))
            # Value is not monotonic at the boundary of the safe zone
            # Hence, this is not an equilibrium
            next_mat = false
            verbose && println()
            verbose && println("    Found no savings equilibrium (value function not monotonic).")
        else
            apLowAtLim = zeros(S, mBase.yParams.yPoints)
            vLowAtLim = zeros(F, mBase.yParams.yPoints)
            apLowAtLimRisky = zeros(S, mBase.yParams.yPoints)
            vLowAtLimRisky = zeros(F, mBase.yParams.yPoints)
            apLowAtLimSafe = zeros(S, mBase.yParams.yPoints)
            vLowAtLimSafe = zeros(F, mBase.yParams.yPoints)

            for yInd in 1:(mBase.yParams.yPoints)
                apLowAtLim[yInd], vLowAtLim[yInd], tempFeasDummy = solveRepayChoice(
                    mMod,
                    s,
                    sExt,
                    yInd,
                    aExtPoints + 1,
                    s.income.mBounds[1],
                    1,
                    mMod.aPoints + aExtPoints
                )
                apLowAtLimRisky[yInd], vLowAtLimRisky[yInd], tempFeasDummy = (
                    solveRepayChoice(
                        mMod,
                        s,
                        sExt,
                        yInd,
                        aExtPoints + 1,
                        s.income.mBounds[1],
                        1,
                        aExtPoints
                    )
                )
                apLowAtLimSafe[yInd], vLowAtLimSafe[yInd], tempFeasDummy = (
                    solveRepayChoice(
                        mMod,
                        s,
                        sExt,
                        yInd,
                        aExtPoints + 1,
                        s.income.mBounds[1],
                        aExtPoints + 1,
                        mMod.aPoints+aExtPoints
                    )
                )
                verbose && println()
            end

            if (minimum(apLowAtLim) >= (aExtPoints+1))
                # policies are pointed in in the Safe Zone 
                # so we have an equilibrium. 
                validSaveEqm[matInd] = true
                next_mat = true
                verbose && println("    Found a savings equilibrium.")
            elseif tempMaxDist < 1e-8
                # policies are pointing out of the Safe Zone, and we have 
                # computed this to a high degree. So this maturity probably 
                # does not have a saving equilibrium.
                next_mat = true
                verbose && println("    Found no savings equilibrium with high precision.")
            else
                # couldn't compute an equilibrium, but the degree of precision 
                # is low. Try with a higher smoothing factor (next wave)
                next_mat = false
                verbose && println("    Found no savings equilibrium, but low precision. Try higher smoothing parameters.")
            end

            maxDevFromSafe = maximum(vLowAtLim .- vLowAtLimSafe)
            maxDevToRisky = maximum(vLowAtLimRisky .- vLowAtLim)
            paramListOut = [
                invMatParamList[matInd],
                lambdaList[matInd],
                modCoupList[matInd],
                aLimit,
                tempMaxDist,
                maxDevFromSafe,
                maxDevToRisky,
                validSaveEqm[matInd]
            ]
            @debug "parameters" invMatParamList[matInd] lambdaList[matInd] modCoupList[matInd] aLimit tempMaxDist maxDevFromSafe maxDevToRisky validSaveEqm[matInd]

            if save_output && validSaveEqm[matInd]
                # we have a valid equilibrium. Let's save it to  disk. 
                writecsv(
                    joinpath(resultDir, "mat_$maturity"*"_apChoicesAtLim.csv"),
                    hcat(apLowAtLimSafe,apLowAtLim,apLowAtLimRisky)
                )
                writecsv(
                    joinpath(resultDir, "mat_$maturity"*"_valuesAtLim.csv"),
                    hcat(vLowAtLimSafe,vLowAtLim,vLowAtLimRisky)
                )
                # println("\n", paramListOut)
                writecsv(
                    joinpath(resultDir, "mat_$maturity"*"_params.csv"),
                    paramListOut
                )
                writeEssentials(s,sExt,"mat_$maturity",resultDir)
            end
        end
        verbose && println()
        if !next_mat 
            if current_wave_ind < total_waves
                current_wave_ind += 1 
                # increase current wave 
                verbose && println("    Trying higher smoothing pars.")
            else 
                verbose && println("    Already at the highest smoothing pars.")
                next_mat = true
            end
        end
        if next_mat 
            (matInd >= matDim) && break
            println("    Moving on to next maturity.")
            # moving on the next maturity index
            matInd = matInd + 1   
            # restarting the wave of smoothing parameters if required
            restart_waves && (current_waves_ind = 1)
        end
    end 
    validEqm = invMatParamList[validSaveEqm]
    notEqm = setdiff(invMatParamList, validEqm)
    return validEqm, notEqm 
end


function find_savings_equilibria(
    matList,
    par,
    yParams, 
    mParams;
    verbose=true,
    save_output=true,
    waves=[
        (mixQFac=0.75, mixVFac=0.75, maxIter=4_000),
        (mixQFac=0.85, mixVFac=0.85, maxIter=15_000),
        (mixQFac=0.98, mixVFac=0.98, maxIter=30_000)
    ],
    restart_waves=true
)
    @unpack (
        aFinePct,
        aLim,
        aPoints,
        aFinePts,
        aExtPoints,
        aExtBoundMag,
        aBoundsSafe,
        resultDirSave
    ) = par

    if verbose
        println("")
        println("Looking for saving equilibria")
    end 

    longBondParams = LongTermBondRDSaveSpec(
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
    iterMatList = matList

    matListSavingsEqm, notEqmMaturities = find_saving_eqm(
        longBondParams,
        iterMatList,
        aFinePct,
        aLim,
        aPoints,
        aFinePts,
        aExtPoints,
        aExtBoundMag,
        par.theta,
        par.hpen0,
        par.hpen1,
        resultDirSave;
        waves=waves,
        restart_waves=restart_waves,
        verbose=verbose,
        save_output=save_output
    )

    if verbose 
        @info "Done. Found a savings eqm for $matListSavingsEqm"
        if !isempty(notEqmMaturities) 
            @warn "But couldn't compute savings eqm for $notEqmMaturities"
        end
    end

    return matListSavingsEqm, notEqmMaturities
end