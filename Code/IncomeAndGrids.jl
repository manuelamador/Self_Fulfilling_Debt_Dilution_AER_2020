#Define versions of the standard normal cdf for use in the tauchen method
std_norm_cdf(x::T) where {T <: Real} = 0.5 * erfc(-x/sqrt(2))

std_norm_cdf(x::Array{T}) where {T <: Real} = 0.5 .* erfc(-x./sqrt(2))


function normCdfDiff(lb::F,ub::F,mu::F,sd::F) where{F<:Real}
    if sign(ub-mu)!=sign(lb-mu)
        return std_norm_cdf((ub-mu)/sd)-std_norm_cdf((lb-mu)/sd)
    elseif (ub-mu)<=zero(F)
        return std_norm_cdf((ub-mu)/sd)-std_norm_cdf((lb-mu)/sd)
    else
        return std_norm_cdf(-(lb-mu)/sd)-std_norm_cdf(-(ub-mu)/sd)
    end
end


function solveMuVec(muTMat::Array{F,2},muSystemSoln::Array{F,1},rho::F) where{F<:Real}

    piStarH=muTMat[2,1]/(muTMat[1,2]+muTMat[2,1])
    piStarL=muTMat[1,2]/(muTMat[1,2]+muTMat[2,1])

    muSystemMat=zeros(2,2)

    #First row says average jump size is muSystemSoln[1]
    muSystemMat[1,1]=-1.0+(piStarH)^(-1)*piStarL*muTMat[2,1]/(1.0-rho*muTMat[2,2])*rho
    #muSystemMat[1,2]=1.0
    muSystemMat[1,2]=(piStarH)^(-1)*piStarL*muTMat[2,1]/(1.0-rho*muTMat[2,2])*(1.0-rho)/muTMat[1,2]

    #Second Row says unconditional mean is muSystemSoln[2]
    muSystemMat[2,1]=piStarL*(1.0+muTMat[2,1]/(1.0-rho*muTMat[2,2])*rho)
    muSystemMat[2,2]=piStarL*muTMat[2,1]/(1.0-rho*muTMat[2,2])*(1.0-rho)/muTMat[1,2]
    #muSystemMat[2,2]=piStarH


    muVec=inv(muSystemMat)*muSystemSoln

    return muVec
end


#The function tauchen(g) generates a discretized version of an AR(1) process using the convention:
#y_t=mu+rho*y_t-1+sqrt(eta2)*e_t
#Its only argument is a variable g of type ar1RDParams.
#Its output is a tuple of 3 objects. They are:
#1. yOut: a one dimensional array of values of the discretized process
#2. yTMatOut: a two dimensional array containing transition probabilities for the process, using the convention that yTMatOut[i,j] is the probability of transitioning to state i conditional on being in state j
#3. yMeanOut: the analytical long run mean of the process

function tauchenRev(g::AR1RDParams{F,S}) where{F<:Real,S<:Integer}
    yARDim=g.yPoints-1
    yStateGrid=zeros(F,g.yPoints)
    yOut=zeros(F,g.yPoints)
    yGrid=zeros(F,yARDim)

    yTMat=zeros(F,g.yPoints,g.yPoints)

    yBar=0.0


    etaSD=sqrt(g.eta2)
    ySD=sqrt(g.eta2/(1.0-g.rho^2))
    yGrid.=collect(LinRange(g.stdSpan[1]*ySD, g.stdSpan[2]*ySD, yARDim))
    yResLow=yGrid[2]-yGrid[1]
    yResHigh=yGrid[end]-yGrid[end-1]

    yDisasterDemeaned=g.muVec[1]-g.muVec[2]


    if g.inflateEndpoints==true
        for i in 1:(yARDim)
            yTMat[1,i]=g.muTMat[2,2]*std_norm_cdf(((yGrid[1] +0.5*yResLow)-g.rho*yGrid[i])/etaSD)
            yTMat[yARDim,i]=g.muTMat[2,2]*(1-std_norm_cdf(((yGrid[yARDim]-0.5*yResHigh)-g.rho*yGrid[i])/etaSD))
            yTMat[yARDim+1,i]=g.muTMat[1,2]

            for j in 2:(yARDim-1)
                yTMat[j,i]=g.muTMat[2,2]*normCdfDiff(0.5*(yGrid[j]+yGrid[j-1]),0.5*(yGrid[j]+yGrid[j+1]),g.rho*yGrid[i],etaSD)
            end
        end
        yTMat[yARDim+1,yARDim+1]=g.muTMat[1,1]
        if g.ergodicReentry==true
            yTMat[1,yARDim+1]=g.muTMat[2,1]*std_norm_cdf((yGrid[1] +0.5*yResLow)/sqrt(g.eta2/(1.0-g.rho^2)))
            yTMat[yARDim,yARDim+1]=g.muTMat[2,1]*(1-std_norm_cdf((yGrid[yARDim]-0.5*yResHigh)/sqrt(g.eta2/(1.0-g.rho^2))))
            for j in 2:(yARDim-1)
                yTMat[j,yARDim+1]=g.muTMat[2,1]*normCdfDiff(0.5*(yGrid[j]+yGrid[j-1]),0.5*(yGrid[j]+yGrid[j+1]),0.0,sqrt(g.eta2/(1.0-g.rho^2)))
            end
        else
            yTMat[1,yARDim+1]=g.muTMat[2,1]*std_norm_cdf(((yGrid[1] +0.5*yResLow)-g.rho*yDisasterDemeaned)/etaSD)
            yTMat[yARDim,yARDim+1]=g.muTMat[2,1]*(1-std_norm_cdf(((yGrid[yARDim]-0.5*yResHigh)-g.rho*yDisasterDemeaned)/etaSD))
            for j in 2:(yARDim-1)
                yTMat[j,yARDim+1]=g.muTMat[2,1]*normCdfDiff(0.5*(yGrid[j]+yGrid[j-1]),0.5*(yGrid[j]+yGrid[j+1]),g.rho*yDisasterDemeaned,etaSD)
            end
        end

    else
        for i in 1:(yARDim)
            yTMat[1,i]=g.muTMat[2,2]*normCdfDiff(yGrid[1]-0.5*yResLow,0.5*(yGrid[1]+yGrid[2]),g.rho*yGrid[i],etaSD)
            yTMat[yARDim,i]=g.muTMat[2,2]*normCdfDiff(0.5*(yGrid[yARDim-1]+yGrid[yARDim]),yGrid[yARDim]+0.5*yResHigh,g.rho*yGrid[i],etaSD)
            yTMat[yARDim+1,i]=g.muTMat[1,2]

            for j in 2:(yARDim-1)
                yTMat[j,i]=g.muTMat[2,2]*normCdfDiff(0.5*(yGrid[j]+yGrid[j-1]),0.5*(yGrid[j]+yGrid[j+1]),g.rho*yGrid[i],etaSD)
            end
        end

        yTMat[yARDim+1,yARDim+1]=g.muTMat[1,1]
        if g.ergodicReentry==true
            yTMat[1,yARDim+1]=g.muTMat[2,1]*normCdfDiff(yGrid[1]-0.5*yResLow,0.5*(yGrid[1]+yGrid[2]),0.0,sqrt(g.eta2/(1.0-g.rho^2)))
            yTMat[yARDim,yARDim+1]=g.muTMat[2,1]*normCdfDiff(0.5*(yGrid[yARDim-1]+yGrid[yARDim]),yGrid[yARDim]+0.5*yResHigh,0.0,sqrt(g.eta2/(1.0-g.rho^2)))
            for j in 2:(yARDim-1)
                yTMat[j,yARDim+1]=g.muTMat[2,1]*normCdfDiff(0.5*(yGrid[j]+yGrid[j-1]),0.5*(yGrid[j]+yGrid[j+1]),0.0,sqrt(g.eta2/(1.0-g.rho^2)))
            end
        else
            yTMat[1,yARDim+1]=g.muTMat[2,1]*normCdfDiff(yGrid[1]-0.5*yResLow,0.5*(yGrid[1]+yGrid[2]),g.rho*yDisasterDemeaned,etaSD)
            yTMat[yARDim,yARDim+1]=g.muTMat[2,1]*normCdfDiff(0.5*(yGrid[yARDim-1]+yGrid[yARDim]),yGrid[yARDim]+0.5*yResHigh,g.rho*yDisasterDemeaned,etaSD)
            for j in 2:(yARDim-1)
                yTMat[j,yARDim+1]=g.muTMat[2,1]*normCdfDiff(0.5*(yGrid[j]+yGrid[j-1]),0.5*(yGrid[j]+yGrid[j+1]),g.rho*yDisasterDemeaned,etaSD)
            end
        end
    end
    yStateGrid.=vcat(yGrid.+g.muVec[2],[g.muVec[1]])
    yOut.=exp.(yStateGrid)

    yTMatSums=sum(yTMat[1:yARDim,:],dims=1)

    for i in 1:(yARDim)
        yTMat[1:yARDim,i].*=yTMatSums[i]^(-1)*g.muTMat[2,2]
    end
    yTMat[1:yARDim,yARDim+1].*=yTMatSums[yARDim+1]^(-1)*g.muTMat[2,1]




    #return yOut,yTMat,yGrid,muIndGrid
    return yOut,yTMat,exp(yBar+0.5*g.eta2/(1.0-g.rho^2))
end

function genStIncDist(
    m::LongTermBondRDSaveSpec,
    income::CombinedIncomeProcessSave,
    tol::Real,
    maxIter::Integer
) 
    oldDist = ones(m.yParams.yPoints) / m.yParams.yPoints
    newDist = zeros(m.yParams.yPoints)

    iCount = 1
    dDist = 1.0 + tol
    while (iCount <= maxIter) && (dDist > tol)
        mul!(newDist, income.yTMat, oldDist)
        dDist = maximum(abs.(newDist - oldDist))
        if (div(iCount, max(div(maxIter, 2), 1)) == 
                (iCount / max(div(maxIter, 2), 1)))
            print(".")
        end
        copyto!(oldDist, newDist)
        iCount += 1
    end
    deflFac = sum(newDist)^(-1)
    return newDist * deflFac
end


function setupYParams(
    mOrig::LongTermBondRDSaveSpec,
    muSystemSoln::Vector,
    tempX::Vector
) 
    piStarH = mOrig.yParams.muTMat[2, 1] / (mOrig.yParams.muTMat[1, 2] +
        mOrig.yParams.muTMat[2, 1])
    piStarL = mOrig.yParams.muTMat[1, 2] / (mOrig.yParams.muTMat[1, 2] + 
        mOrig.yParams.muTMat[2, 1])

    muSystemMat = zeros(2,2)

    #First row says average jump size is muSystemSoln[1]
    muSystemMat[1, 1] = (
        -1.0 + (piStarH)^(-1) * piStarL * mOrig.yParams.muTMat[2, 1]/
        (1.0 - tempX[1] * mOrig.yParams.muTMat[2, 2]) * tempX[1]
    )
    muSystemMat[1, 2] = (
        (piStarH)^(-1) * piStarL * mOrig.yParams.muTMat[2, 1] /
        (1.0 - tempX[1] * mOrig.yParams.muTMat[2, 2]) * (1.0 - tempX[1]) / 
        mOrig.yParams.muTMat[1, 2]
    )

    #Second Row says unconditional mean is muSystemSoln[2]
    muSystemMat[2, 1] = (
        piStarL * (1.0 + mOrig.yParams.muTMat[2, 1] / 
        (1.0 - tempX[1] * mOrig.yParams.muTMat[2, 2]) * tempX[1])
    )
    muSystemMat[2, 2] = (
        piStarL * mOrig.yParams.muTMat[2, 1] /
        (1.0 - tempX[1] * mOrig.yParams.muTMat[2, 2]) * (1.0 - tempX[1]) / 
        mOrig.yParams.muTMat[1, 2]
    )

    muVec = inv(muSystemMat)*muSystemSoln

    return AR1RDParams(
        mOrig.yParams.yPoints,
        tempX[1],
        tempX[2]^2,
        muVec,
        mOrig.yParams.muTMat,
        mOrig.yParams.stdSpan,
        mOrig.yParams.inflateEndpoints,
        mOrig.yParams.ergodicReentry
    )
end


function makeEmpiricalMoments(
    mNew::LongTermBondRDSaveSpec,
    income::CombinedIncomeProcessSave,
    targetRho::Real,
    targetEta::Real;
    checkVar::Bool=false
) 
    stDistY = genStIncDist(mNew, income, 1e-30, 10000)
    logY = log.(income.yGrid)
    logYMean = dot(stDistY,logY)
    varY = dot(stDistY,(logY.-logYMean).^2)
    yAutoCov = (
        dot(stDistY,sum(income.yTMat.*(logY * logY'), dims=1)) - logYMean^2
    )
    rhoYMoment = yAutoCov / varY
    eta2YMoment = dot(
        stDistY, 
        sum(income.yTMat.*(logY.-rhoYMoment.*logY').^2, dims=1)
    )
    etaYMoment = sqrt(eta2YMoment)
    if checkVar == false
        return rhoYMoment - targetRho, etaYMoment - targetEta
    else
        return (
            rhoYMoment - targetRho, 
            etaYMoment - targetEta, 
            sqrt(varY) - sqrt(targetEta^2 / (1.0 - targetRho^2))
        )
    end
end


function getStepDirection(
    tempF::Vector,
    tempFAtLB::Vector,
    tempFAtUB::Vector,
    tempLBX::Vector,
    tempUBX::Vector
) 
    stepInRhoSpace = false
    stepInEtaSpace = false

    tempXDiff = tempUBX.-tempLBX

    alwaysStepInRhoSpace = ifelse(
        (tempXDiff[1] > 10 * tempXDiff[2]) && (tempF[1] !=0) && (
            (sign(tempFAtLB[1]) != sign(tempF[1])) || 
            (sign(tempFAtUB[1]) != sign(tempF[1]))
            ),
        true,
        false
    )
    alwaysStepInEtaSpace = ifelse(
        (tempXDiff[2] > 10 * tempXDiff[1]) && (tempF[2]!=0) && (
            (sign(tempFAtLB[2]) != sign(tempF[2])) || 
            (sign(tempFAtUB[2]) != sign(tempF[2]))
            ),
        true,
        false
    )

    if (
        (abs(tempF[1]) > abs(tempF[2])) && 
        (alwaysStepInRhoSpace == false) && 
        (alwaysStepInEtaSpace == false)
    )
        if (
            (sign(tempFAtLB[1]) != sign(tempF[1])) ||
            (sign(tempFAtUB[1]) != sign(tempF[1]))
        )
            stepInRhoSpace = true
        elseif (
            (sign(tempFAtLB[2]) != sign(tempF[2])) ||
            (sign(tempFAtUB[2]) != sign(tempF[2]))
        )
            stepInEtaSpace=true
        end
    elseif (
        (abs(tempF[2]) != 0.0) &&
        (alwaysStepInRhoSpace == false) &&
        (alwaysStepInEtaSpace == false)
    )
        if (
            (sign(tempFAtLB[2]) != sign(tempF[2])) ||
            (sign(tempFAtUB[2]) != sign(tempF[2]))
        )
            stepInEtaSpace = true
        elseif (
            (sign(tempFAtLB[1]) != sign(tempF[1])) ||
            (sign(tempFAtUB[1])!=sign(tempF[1]))
        )
            stepInRhoSpace = true
        end
    elseif alwaysStepInRhoSpace == true
        if (
            (sign(tempFAtLB[1]) != sign(tempF[1])) ||
            (sign(tempFAtUB[1])!=sign(tempF[1]))
        )
            stepInRhoSpace = true
        end
    elseif alwaysStepInEtaSpace == true
        if (
            (sign(tempFAtLB[2]) != sign(tempF[2])) ||
            (sign(tempFAtUB[2])!=sign(tempF[2]))
        )
            stepInEtaSpace=true
        end
    end
    stepTowardsLB = false
    stepTowardsUB = false
    if stepInRhoSpace == true
        if tempFAtLB[1] < 0.0
            if tempF[1] > 0.0
                stepTowardsLB = true
            else
                stepTowardsUB = true
            end
        else
            if tempF[1]>0.0
                stepTowardsUB = true
            else
                stepTowardsLB = true
            end
        end
    elseif stepInEtaSpace == true
        if tempFAtLB[2] < 0.0
            if tempF[2] > 0.0
                stepTowardsLB = true
            else
                stepTowardsUB = true
            end
        else
            if tempF[2] > 0.0
                stepTowardsUB = true
            else
                stepTowardsLB = true
            end
        end

    end

    return stepInRhoSpace, stepInEtaSpace, stepTowardsLB, stepTowardsUB
end


function fitRDMomentsAR1(
    mOrig::LongTermBondRDSaveSpec{F, S},
    RDMag::F,
    targetRho::F,
    targetEta::F,
    rhoBounds::Array{F, 1},
    etaBounds::Array{F, 1},
    tol::F,
    maxIter::S
) where{F<:Real, S<:Integer}

    tempX = zeros(F,2)      
    tempF = zeros(F,2)

    muSystemSoln = [RDMag, 0.0]

    tempFAtLB = zeros(F, 2)
    tempFAtUB = zeros(F, 2)

    tempLBX = [rhoBounds[1], etaBounds[1]]
    tempUBX = [rhoBounds[2], etaBounds[2]]
    tempX .= tempLBX

    yParamsNew = setupYParams(mOrig, muSystemSoln, tempX)
    mNew = LongTermBondRDSaveSpec(
        mOrig.beta,
        mOrig.gamma,
        mOrig.R,
        mOrig.lambda,
        mOrig.coup,
        mOrig.aBounds,
        mOrig.aPoints,
        yParamsNew,
        mOrig.mParams
    )
    income = makeIncomeProcess(mNew)

    tempFAtLB .= makeEmpiricalMoments(mNew, income, targetRho, targetEta)

    tempX .= tempUBX
    yParamsNew = setupYParams(mOrig, muSystemSoln, tempX)
    mNew = LongTermBondRDSaveSpec(
        mOrig.beta,
        mOrig.gamma,
        mOrig.R,
        mOrig.lambda,
        mOrig.coup,
        mOrig.aBounds,
        mOrig.aPoints,
        yParamsNew,
        mOrig.mParams
    )
    income = makeIncomeProcess(mNew)

    tempFAtUB .= makeEmpiricalMoments(mNew, income, targetRho, targetEta)
    tempX = 0.5 .* (tempLBX .+ tempUBX)
    yParamsNew = setupYParams(mOrig, muSystemSoln, tempX)
    mNew = LongTermBondRDSaveSpec(
        mOrig.beta,
        mOrig.gamma,
        mOrig.R,
        mOrig.lambda,
        mOrig.coup,
        mOrig.aBounds,
        mOrig.aPoints,
        yParamsNew,
        mOrig.mParams
    )
    income = makeIncomeProcess(mNew)

    tempF .= makeEmpiricalMoments(mNew, income, targetRho, targetEta)

    (
        stepInRhoSpace,
        stepInEtaSpace,
        stepTowardsLB,
        stepTowardsUB
    ) = getStepDirection(tempF, tempFAtLB, tempFAtUB, tempLBX, tempUBX)

    iCount = 2

    bestSolnYet = deepcopy(tempX)
    bestMomentsYet = deepcopy(tempF)

    if maximum(abs.(tempFAtLB)) < maximum(abs.(bestMomentsYet))
        bestSolnYet .= tempLBX
        bestMomentsYet .= tempFAtLB
    end
    if maximum(abs.(tempFAtUB)) < maximum(abs.(bestMomentsYet))
        bestSolnYet .= tempUBX
        bestMomentsYet .= tempFAtUB
    end

    # println([iCount-1,tempLBX[1],tempFAtLB[1],tempX[1],tempF[1],tempF[1]+targetRho,tempUBX[1],tempFAtUB[1]])
    # println([iCount-1,tempLBX[2],tempFAtLB[2],tempX[2],tempF[2],tempF[2]+targetEta,tempUBX[2],tempFAtUB[2]])

    while (
        (iCount <= maxIter) && 
        (maximum(abs.(tempF)) > tol) &&
        ((stepInRhoSpace == true) || (stepInEtaSpace == true)) &&
        (maximum(abs.(bestMomentsYet)) > tol)
    )
        if stepInRhoSpace == true
            if stepTowardsLB == true
                tempUBX[1] = 0.5 * (tempUBX[1] + tempX[1])
            else
                tempLBX[1] = 0.5 * (tempLBX[1] + tempX[1])
            end
            tempX[1] = 0.5 * (tempLBX[1] + tempUBX[1])
        elseif stepInEtaSpace == true
            if stepTowardsLB == true
                tempUBX[2] = 0.5 * (tempUBX[2] + tempX[2])
            else
                tempLBX[2] = 0.5 * (tempLBX[2] + tempX[2])
            end
            tempX[2] = 0.5 * (tempLBX[2] + tempUBX[2])
        end

        if stepTowardsLB == true
            yParamsNew = setupYParams(mOrig, muSystemSoln, tempUBX)
            mNew = LongTermBondRDSaveSpec(
                mOrig.beta,
                mOrig.gamma,
                mOrig.R,
                mOrig.lambda,
                mOrig.coup,
                mOrig.aBounds,
                mOrig.aPoints,
                yParamsNew,
                mOrig.mParams
            )
            income=makeIncomeProcess(mNew)

            tempFAtUB .= makeEmpiricalMoments(
                mNew, income, targetRho, targetEta
            )
            if maximum(abs.(tempFAtUB)) < maximum(abs.(bestMomentsYet))
                bestSolnYet .= tempUBX
                bestMomentsYet .= tempFAtUB
            end
        elseif stepTowardsUB == true
            yParamsNew = setupYParams(mOrig, muSystemSoln, tempLBX)
            mNew=LongTermBondRDSaveSpec(
                mOrig.beta,
                mOrig.gamma,
                mOrig.R,
                mOrig.lambda,
                mOrig.coup,
                mOrig.aBounds,
                mOrig.aPoints,
                yParamsNew,
                mOrig.mParams
            )
            income = makeIncomeProcess(mNew)
            tempFAtLB .= makeEmpiricalMoments(mNew, income, targetRho, targetEta)
            if maximum(abs.(tempFAtLB)) < maximum(abs.(bestMomentsYet))
                bestSolnYet .= tempLBX
                bestMomentsYet .= tempFAtLB
            end
        end

        yParamsNew = setupYParams(mOrig, muSystemSoln, tempX)
        mNew = LongTermBondRDSaveSpec(
            mOrig.beta,
            mOrig.gamma,
            mOrig.R,
            mOrig.lambda,
            mOrig.coup,
            mOrig.aBounds,
            mOrig.aPoints,
            yParamsNew,
            mOrig.mParams
        )
        income = makeIncomeProcess(mNew)

        tempF .= makeEmpiricalMoments(mNew, income, targetRho, targetEta)
        if maximum(abs.(tempF)) < maximum(abs.(bestMomentsYet))
            bestSolnYet .= tempX
            bestMomentsYet .= tempF
        end
        (
            stepInRhoSpace,
            stepInEtaSpace,
            stepTowardsLB,
            stepTowardsUB
        ) = getStepDirection(tempF, tempFAtLB, tempFAtUB, tempLBX, tempUBX)

        # println([iCount,tempLBX[1],tempFAtLB[1],tempX[1],tempF[1],tempF[1]+targetRho,tempUBX[1],tempFAtUB[1]])
        # println([iCount,tempLBX[2],tempFAtLB[2],tempX[2],tempF[2],tempF[2]+targetEta,tempUBX[2],tempFAtUB[2]])
        print(".")
        
        iCount += 1
    end

    yParamsNew = setupYParams(mOrig, muSystemSoln, bestSolnYet)
    mNew = LongTermBondRDSaveSpec(
        mOrig.beta,
        mOrig.gamma,
        mOrig.R,
        mOrig.lambda,
        mOrig.coup,
        mOrig.aBounds,
        mOrig.aPoints,
        yParamsNew,
        mOrig.mParams
    )
    income = makeIncomeProcess(mNew)
    momentsOut = zeros(F,3)
    momentsOut .= makeEmpiricalMoments(
        mNew,
        income,
        targetRho,
        targetEta,
        checkVar=true
    )

    return (rho=bestSolnYet[1], eta=bestSolnYet[2]), momentsOut
end


function obtain_income_process(par; verbose=true)
    @unpack (
        beta,
        theta,
        gamma,
        R,
        delta,
        kappa,
        y_points,
        rho_target,
        eta_squared_target,
        std_span,
        inflate_end_points,
        ergodic_rentry,
        mu_T,
        mu_system_soln,
        m_points,
        eps_squared,
        m_mu,
        m_std_span,
    ) = par 

    if verbose
        println("Finding the parameters for AR1")
    end 
    
    coup = kappa
    lambda = delta

    mu_vec = solveMuVec(mu_T, mu_system_soln, rho_target)
    y_params = AR1RDParams(
        y_points,
        rho_target,
        eta_squared_target,
        mu_vec,
        mu_T,
        std_span,
        inflate_end_points,
        ergodic_rentry
    )
    
    m_params = IIDParams(
        m_points,
        eps_squared,
        m_mu,
        m_std_span
    )
    
    long_bond_params = LongTermBondRDSaveSpec(
        beta,
        gamma,
        R,
        lambda,
        coup,
        [-1.0, 0.0], # abounds  (not needed for this)
        1,  # apoints (not needed for this)
        y_params,
        m_params
    )
    # ^ Constructs a struct with the model parameters
        
    fitted_rho_eta, = fitRDMomentsAR1(
        long_bond_params,
        mu_system_soln[1],
        rho_target, # target values
        sqrt(eta_squared_target), # target_values
        [0.9, 0.98],  # rho bounds
        [0.015, 0.035],  # eta bounds
        1e-16, # tolerance
        1000 # maxIters
    )
    
    # update rho to its fitted value
    rho_fitted = fitted_rho_eta.rho
    eta_squared_fitted = fitted_rho_eta.eta^2
      
    # generate the new income process parameters with the fitted 
    # values of rho and eta
    mu_vec = solveMuVec(mu_T, mu_system_soln, rho_fitted)

    y_params = AR1RDParams(
        y_points,
        rho_fitted,
        eta_squared_fitted,
        mu_vec,
        mu_T,
        std_span,
        inflate_end_points,
        ergodic_rentry
    )

    if verbose 
        println("")
        @info "Fitted values of rho and eta $fitted_rho_eta"
    end

    return y_params, m_params, fitted_rho_eta
end




#genBisectSearchGrids is a utility function for constructing a gridSearchParams object
#Its only argument is the number of points in the grid of interest
#Its output is simply an object of type gridSearchParams
function genBisectSearchGrids(numPoints::S) where{S<:Integer}
    #Initialize all the components of the output object
    numLevels=S(ceil(log(numPoints-1)/log(2)))
    levelPoints=Array{S,1}(undef,numLevels)
    boundGridPoints=Array{S,1}(undef,numLevels)

    boundGrids=Array{Array{S,1},1}(undef,numLevels)
    pointGrids=Array{Array{S,1},1}(undef,numLevels)
    pointLB=Array{Array{S,1},1}(undef,numLevels)
    pointUB=Array{Array{S,1},1}(undef,numLevels)

    #Initialize a uncleaned version of boundGrids
    boundGridsFull=Array{Array{S,1},1}(undef,numLevels)
    boundGridsFull[1]=[1,numPoints]

    #Set the values of each objects at the first level
    boundGrids[1]=[1,numPoints]
    pointGrids[1]=[S(ceil(numPoints/2))]
    levelPoints[1]=1
    boundGridPoints[1]=2
    pointLB[1]=[1]
    pointUB[1]=[numPoints]

    #Iterate over levels
    for j in 2:numLevels
        #Set the current number of points to 0
        levelPoints[j]=0

        #Construct the set of points available to be used as upper or lower bounds
        tempBoundGrid=vcat(boundGridsFull[j-1],pointGrids[j-1])

        #Sort this array and remove all duplicates
        sort!(tempBoundGrid,alg=QuickSort)
        notFinishedCleaning=true
        tempInd=2
        while notFinishedCleaning==true
            if tempBoundGrid[tempInd]==tempBoundGrid[tempInd-1]
                deleteat!(tempBoundGrid,tempInd)
                tempInd==2
            else
                tempInd+=1
            end
            if tempInd>length(tempBoundGrid)
                notFinishedCleaning=false
            end
        end
        #Record the current list of available bounds
        boundGridsFull[j]=deepcopy(tempBoundGrid)
        #Initialize a set of arrays which will mark whether or not the point in the temporary bound grid is actually needed to calculate any of the new points
        tempKeepMark=ones(Bool,length(tempBoundGrid))
        tempKeepMark[1]=ifelse((tempBoundGrid[1]+1)==tempBoundGrid[2],false,true)
        tempKeepMark[end]=ifelse((tempBoundGrid[end-1]+1)==tempBoundGrid[end],false,true)
        #If necessary, iterate over the elements of the temporary bound grid
        #Check whether each interior element is both 1 greater than the previous element and 1 less than the next element. If this is true, the point is not a bound index for any points at this level
        #Also check whether each element except the first is equal to the element preceding it. Whenever this is the case, a new point will be generated, and that must be reflected by levelPoints

        if length(tempBoundGrid)>2
            for i in 2:(length(tempBoundGrid)-1)
                if ((tempBoundGrid[i-1])==(tempBoundGrid[i]-1))&&((tempBoundGrid[i+1])==(tempBoundGrid[i]+1))
                    tempKeepMark[i]=false
                end
            end
            for i in 2:(length(tempBoundGrid))
                if ((tempBoundGrid[i-1])!=(tempBoundGrid[i]-1))
                    levelPoints[j]+=1
                end
            end
        end
        #Generate the the current boundGrid, record its length, and initialize the grid of point indices and lower and upper bound indices at this level
        boundGrids[j]=tempBoundGrid[tempKeepMark]
        boundGridPoints[j]=length(boundGrids[j])
        pointGrids[j]=zeros(S,levelPoints[j])
        pointLB[j]=zeros(S,levelPoints[j])
        pointUB[j]=zeros(S,levelPoints[j])

        #Initialize a variable tracking the index in pointGrids[j] next point to be added
        tempInd=1
        #Alternate between using ceilings and floors
        if isodd(j)==true
            #Iterate over the elements of the temporary bound grid
            for k in 1:(length(tempBoundGrid)-1)
                #Check whether a point needs to be added and, if necessary, add it and record the relevant upper and lower bound indices
                if tempBoundGrid[k]!=(tempBoundGrid[k+1]-1)
                    pointGrids[j][tempInd]=S(ceil((tempBoundGrid[k]+tempBoundGrid[k+1])/2))
                    pointLB[j][tempInd]=tempBoundGrid[k]
                    pointUB[j][tempInd]=tempBoundGrid[k+1]

                    tempInd+=1
                end
            end
        else
            #Iterate over the elements of the temporary bound grid
            for k in 1:(length(tempBoundGrid)-1)
                #Check whether a point needs to be added and, if necessary, add it and record the relevant upper and lower bound indices
                if tempBoundGrid[k]!=(tempBoundGrid[k+1]-1)
                    pointGrids[j][tempInd]=S(floor((tempBoundGrid[k]+tempBoundGrid[k+1])/2))
                    pointLB[j][tempInd]=tempBoundGrid[k]
                    pointUB[j][tempInd]=tempBoundGrid[k+1]

                    tempInd+=1
                end
            end
        end
    end
    return GridSearchParams(boundGrids,pointGrids,pointLB,pointUB,boundGridPoints,levelPoints,numLevels)
end

