function makeConsPolicy(m::LongTermBondRDBorrowSpec{F,S}) where{F<:Real,S<:Integer}
    consRGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    consGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    consNoPen=zeros(F,m.aPoints,m.yParams.yPoints)
    consNoPenD=zeros(F,m.yParams.yPoints)
    consRVar=zeros(F,m.aPoints,m.yParams.yPoints)
    consVar=zeros(F,m.aPoints,m.yParams.yPoints)
    consDVar=zeros(F,m.yParams.yPoints)

    return ConsumptionPolicy(consRGrid,consGrid,consNoPen,consNoPenD,consRVar,consVar,consDVar)

end


function makeConsDivYPolicy(m::LongTermBondRDBorrowSpec{F,S}) where{F<:Real,S<:Integer}
    consDivYGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    consRDivYGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    consDivYDInitGrid=zeros(F,m.yParams.yPoints)
    consDivYDFutGrid=zeros(F,m.yParams.yPoints)
    consRDivYVar=zeros(F,m.aPoints,m.yParams.yPoints)
    consDivYVar=zeros(F,m.aPoints,m.yParams.yPoints)
    consDivYDFutVar=zeros(F,m.yParams.yPoints)
    consDivYDInitVar=zeros(F,m.yParams.yPoints)

    return ConsumptionDivYPolicy(consDivYGrid,consRDivYGrid,consDivYDInitGrid,consDivYDFutGrid,consRDivYVar,consDivYVar,consDivYDFutVar,consDivYDInitVar)

end
function makeUtilityPolicy(m::LongTermBondRDBorrowSpec{F,S}) where{F<:Real,S<:Integer}
    uGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    uRGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    uDInitGrid=zeros(F,m.yParams.yPoints)
    uDFutGrid=zeros(F,m.yParams.yPoints)
    uNoPen=zeros(F,m.aPoints,m.yParams.yPoints)
    uNoPenD=zeros(F,m.yParams.yPoints)
    utYObs=zeros(F,m.aPoints,m.yParams.yPoints)
    utYObsD=zeros(F,m.yParams.yPoints)

    return UtilityPolicy(uGrid,uRGrid,uDInitGrid,uDFutGrid,uNoPen,uNoPenD,utYObs,utYObsD)
end

function makeOutputPolicy(m::LongTermBondRDBorrowSpec{F,S}) where{F<:Real,S<:Integer}
    yObs=zeros(F,m.aPoints,m.yParams.yPoints)
    yObsD=zeros(F,m.yParams.yPoints)

    return OutputPolicy(yObs,yObsD)
end

function makeAPDefPolicy(m::LongTermBondRDBorrowSpec{F,S}) where{F<:Real,S<:Integer}
    apRGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    apRDivYGrid=zeros(F,m.aPoints,m.yParams.yPoints)
    defaultProb=zeros(F,m.aPoints,m.yParams.yPoints)
    repayProb=zeros(F,m.aPoints,m.yParams.yPoints)
    apProbability=Array{SparseMatrixCSC{F,S},1}(undef,m.yParams.yPoints)
    lastAlwaysDefInd=zeros(S,m.yParams.yPoints)

    SDummyF=sparse(zeros(F,m.aPoints,m.aPoints))
    for i in 1:(m.yParams.yPoints)
        apProbability[i]=deepcopy(SDummyF)
    end

    return APrimeDefPolicy(apRGrid,apRDivYGrid,defaultProb,repayProb,apProbability,lastAlwaysDefInd)

end


function makeMeanPolicies(m::LongTermBondRDBorrowSpec{F,S},s::LongTermBondRDBorrowEval{F,S,T}) where{F<:Real,S<:Integer,T}


    cPol=makeConsPolicy(m)
    cDivYPol=makeConsDivYPolicy(m)
    apDefPol=makeAPDefPolicy(m)
    uPol=makeUtilityPolicy(m)
    yPol=makeOutputPolicy(m)

    copyto!(uPol.uDInitGrid,s.VF.vDInitFlow)
    copyto!(uPol.uDFutGrid,s.VF.vDFutFlow)

    copyto!(cPol.consNoPenD,s.income.yGrid)
    copyto!(yPol.yObsD,s.income.yDefGrid)
    copyto!(uPol.utYObsD,s.VF.vDFutFlow)


    for i in 1:m.yParams.yPoints
        tempADInd=0
        while tempADInd<m.aPoints
            if s.pol.alwaysDefault[tempADInd+1,i]==false
                break
            else
                tempADInd+=1
            end
        end
        apDefPol.lastAlwaysDefInd[i]=tempADInd
    end


    for i in 1:(m.yParams.yPoints-1)
        tempPSum=0.0
        for mInd in 1:(m.mParams.mPoints-1)
            cDivYPol.consDivYDInitGrid[i]+=s.income.mProb[mInd]*(s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])
            cDivYPol.consDivYDFutGrid[i]+=s.income.mProb[mInd]*(s.income.yDefGrid[i]+s.income.mMidPoints[mInd])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])
            uPol.uNoPenD[i]+=s.income.mProb[mInd]*u(m,s.income.yGrid[i]+s.income.mMidPoints[mInd])
            tempPSum+=s.income.mProb[mInd]
        end
        tempPSumInv=tempPSum^(-1)
        cDivYPol.consDivYDInitGrid[i]*=tempPSumInv
        cDivYPol.consDivYDFutGrid[i]*=tempPSumInv
        uPol.uNoPenD[i]*=tempPSumInv
    end
    cDivYPol.consDivYDInitGrid[end]=s.income.yDefGrid[end]/s.income.yDefGrid[end]
    cDivYPol.consDivYDFutGrid[end]=s.income.yDefGrid[end]/s.income.yDefGrid[end]

    uPol.uNoPenD[end]=u(m,s.income.yGrid[end])

    for i in 1:(m.yParams.yPoints-1)
        tempPSum=0.0
        for mInd in 1:(m.mParams.mPoints-1)
            cPol.consDVar[i]+=s.income.mProb[mInd]*(s.income.yDefGrid[i]+s.income.mMidPoints[mInd]-s.income.yDefGrid[i])^2
            cDivYPol.consDivYDFutVar[i]+=s.income.mProb[mInd]*((s.income.yDefGrid[i]+s.income.mMidPoints[mInd])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])-cDivYPol.consDivYDFutGrid[i])^2
            cDivYPol.consDivYDInitVar[i]+=s.income.mProb[mInd]*((s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mInd])-cDivYPol.consDivYDInitGrid[i])^2
            tempPSum+=s.income.mProb[mInd]
        end
        tempPSumInv=tempPSum^(-1)
        cPol.consDVar[i]*=tempPSumInv
        cDivYPol.consDivYDFutVar[i]*=tempPSumInv
        cDivYPol.consDivYDInitVar[i]*=tempPSumInv

    end
    cPol.consDVar[end]=0.0
    cDivYPol.consDivYDFutVar[end]=1.0
    cDivYPol.consDivYDInitVar[end]=1.0


    for i in 1:(m.yParams.yPoints-1)
        for j in 1:(m.aPoints)
            tempPSum=0.0
            if s.pol.neverDefault[j,i]==true
                mUBInd=2
                aPolicyInd=1
                lastMVal=s.income.mGrid[1]
                while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                tempPSumInv=tempPSum^(-1)
                cPol.consRGrid[j,i]*=tempPSumInv
                apDefPol.apRGrid[j,i]*=tempPSumInv
                apDefPol.apRDivYGrid[j,i]*=tempPSumInv
                cDivYPol.consRDivYGrid[j,i]*=tempPSumInv
                uPol.uRGrid[j,i]*=tempPSumInv
                for apInd in 1:(s.pol.mListLength[j,i])
                    apDefPol.apProbability[i][apInd,j]*=tempPSumInv
                end
                cPol.consGrid[j,i]=cPol.consRGrid[j,i]
                cDivYPol.consDivYGrid[j,i]=cDivYPol.consRDivYGrid[j,i]
                uPol.uGrid[j,i]=uPol.uRGrid[j,i]

                uPol.utYObs[j,i]=uPol.uNoPenD[i]
                cPol.consNoPen[j,i]=cPol.consGrid[j,i]
                uPol.uNoPen[j,i]=uPol.uGrid[j,i]
                yPol.yObs[j,i]=s.income.yGrid[i]
            elseif s.pol.alwaysDefault[j,i]==false
                #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
                #and set the last value of m observed to the value directly preceding that upper bound
                mUBInd=searchsortedfirst(s.income.mGrid,s.pol.defThreshold[j,i])
                lastMVal=s.income.mGrid[mUBInd-1]
                #If there are any entire intervals in which default occurs, add their contribution to the output value
                if mUBInd>2
                    for k in 3:mUBInd
                        apDefPol.defaultProb[j,i]+=s.income.mProb[k-2]
                        tempPSum+=s.income.mProb[k-2]
                        cDivYPol.consDivYGrid[j,i]+=s.income.mProb[k-2]*(s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[k-2])

                        uPol.utYObs[j,i]+=s.income.mProb[k-2]*u(m,s.income.yDefGrid[i]+s.income.mBounds[1])
                        yPol.yObs[j,i]+=s.income.mProb[k-2]*(s.income.yDefGrid[i]+s.income.mBounds[1])
                        uPol.uNoPen[j,i]+=s.income.mProb[k-2]*u(m,s.income.yGrid[i]+s.income.mMidPoints[k-2])
                        cPol.consNoPen[j,i]+=s.income.mProb[k-2]*(s.income.yGrid[i]+s.income.mMidPoints[k-2])
                    end
                end
                #Set the index of the first relevant entry of the borowing policy function
                aPolicyInd=s.pol.firstRepayInd[j,i]

                #Add the contribution of default in the interval in which the threshold lies to the output value
                apDefPol.defaultProb[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes
                cDivYPol.consDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])


                uPol.utYObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*u(m,s.income.yDefGrid[i]+s.income.mBounds[1])
                yPol.yObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yDefGrid[i]+s.income.mBounds[1])
                uPol.uNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*u(m,s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                cPol.consNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])

                tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes

                #Update the last value of m observed to the threshold level of m at which default occurs
                lastMVal=s.pol.defThreshold[j,i]

                while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])

                        uPol.utYObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        yPol.yObs[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*u(m,s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])
                        cPol.consNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        apDefPol.apProbability[i][aPolicyInd,j]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        cPol.consRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])
                        apDefPol.apRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]
                        apDefPol.apRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*s.aGrid[s.pol.apPolicy[i][aPolicyInd,j]]/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        cDivYPol.consRDivYGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uRGrid[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])

                        uPol.utYObs[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        yPol.yObs[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])
                        uPol.uNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*u(m,s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])
                        cPol.consNoPen[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])


                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                cDivYPol.consDivYGrid[j,i]+=cDivYPol.consRDivYGrid[j,i]

                tempPSumInv=tempPSum^(-1)
                apDefPol.defaultProb[j,i]*=tempPSumInv
                cDivYPol.consDivYGrid[j,i]*=tempPSumInv

                uPol.utYObs[j,i]*=tempPSumInv
                yPol.yObs[j,i]*=tempPSumInv
                uPol.uNoPen[j,i]*=tempPSumInv
                cPol.consNoPen[j,i]*=tempPSumInv

                cPol.consRGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                apDefPol.apRGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                apDefPol.apRDivYGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                cDivYPol.consRDivYGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                uPol.uRGrid[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)

                cPol.consGrid[j,i]=(1.0-apDefPol.defaultProb[j,i])*cPol.consRGrid[j,i]+apDefPol.defaultProb[j,i]*(s.income.yDefGrid[i]+s.income.mBounds[1])
                uPol.uGrid[j,i]=(1.0-apDefPol.defaultProb[j,i])*uPol.uRGrid[j,i]+apDefPol.defaultProb[j,i]*uPol.uDInitGrid[i]

                for apInd in (s.pol.firstRepayInd[j,i]):(s.pol.mListLength[j,i])
                    apDefPol.apProbability[i][apInd,j]*=tempPSumInv
                end
            else
                apDefPol.defaultProb[j,i]=1.0
                cPol.consGrid[j,i]=(s.income.yDefGrid[i]+s.income.mBounds[1])
                cDivYPol.consDivYGrid[j,i]=cDivYPol.consDivYDInitGrid[i]
                uPol.uGrid[j,i]=uPol.uDInitGrid[i]
                uPol.utYObs[j,i]=u(m,s.income.yDefGrid[i]+s.income.mBounds[1])
                yPol.yObs[j,i]=s.income.yDefGrid[i]+s.income.mBounds[1]
                uPol.uNoPen[j,i]=uPol.uNoPenD[i]
                cPol.consNoPen[j,i]=cPol.consNoPenD[i]
            end

        end
    end
    for j in 1:m.aPoints
        if s.pol.alwaysDefault[j,end]==false

            apDefPol.apProbability[end][1,j]=1.0
            cPol.consRGrid[j,end]=s.consM0Grid[s.pol.apLowM[j,end],j,end]
            apDefPol.apRGrid[j,end]=s.aGrid[s.pol.apLowM[j,end]]
            apDefPol.apRDivYGrid[j,end]=s.aGrid[s.pol.apLowM[j,end]]/s.income.yGrid[end]
            cDivYPol.consRDivYGrid[j,end]=s.consM0Grid[s.pol.apLowM[j,end],j,end]/s.income.yGrid[end]
            uPol.uRGrid[j,end]=u(m,s.consM0Grid[s.pol.apLowM[j,end],j,end])


            cPol.consGrid[j,end]=cPol.consRGrid[j,end]
            cDivYPol.consDivYGrid[j,end]=cDivYPol.consRDivYGrid[j,end]
            uPol.uGrid[j,end]=uPol.uRGrid[j,end]

            uPol.utYObs[j,end]=uPol.uNoPenD[end]
            cPol.consNoPen[j,end]=cPol.consGrid[j,end]
            uPol.uNoPen[j,end]=uPol.uGrid[j,end]
            yPol.yObs[j,end]=s.income.yGrid[end]

            cPol.consRVar[j,end]=0.0
            cDivYPol.consRDivYVar[j,end]=0.0
            cPol.consVar[j,end]=0.0
            cDivYPol.consDivYVar[j,end]=0.0
        else
            apDefPol.defaultProb[j,end]=1.0
            cPol.consGrid[j,end]=s.income.yDefGrid[end]
            cDivYPol.consDivYGrid[j,end]=cDivYPol.consDivYDInitGrid[end]
            uPol.uGrid[j,end]=uPol.uDInitGrid[end]
            uPol.utYObs[j,end]=u(m,s.income.yDefGrid[end])
            yPol.yObs[j,end]=s.income.yDefGrid[end]
            uPol.uNoPen[j,end]=uPol.uNoPenD[end]
            cPol.consNoPen[j,end]=cPol.consNoPenD[end]
            cPol.consRVar[j,end]=0.0
            cDivYPol.consRDivYVar[j,end]=0.0
            cPol.consVar[j,end]=0.0
            cDivYPol.consDivYVar[j,end]=0.0
        end
    end


    for i in 1:(m.yParams.yPoints-1)
        for j in 1:(m.aPoints)
            tempPSum=0.0
            if s.pol.neverDefault[j,i]==true
                mUBInd=2
                aPolicyInd=1
                lastMVal=s.income.mGrid[1]
                while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1]-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*((s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1]-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*((s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                tempPSumInv=tempPSum^(-1)
                cPol.consRVar[j,i]*=tempPSumInv
                cDivYPol.consRDivYVar[j,i]*=tempPSumInv

                cPol.consVar[j,i]=cPol.consRVar[j,i]
                cDivYPol.consDivYVar[j,i]=cDivYPol.consRDivYVar[j,i]

            elseif s.pol.alwaysDefault[j,i]==false
                #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
                #and set the last value of m observed to the value directly preceding that upper bound
                mUBInd=searchsortedfirst(s.income.mGrid,s.pol.defThreshold[j,i])
                lastMVal=s.income.mGrid[mUBInd-1]
                #If there are any entire intervals in which default occurs, add their contribution to the output value
                if mUBInd>2
                    for k in 3:mUBInd
                        tempPSum+=s.income.mProb[k-2]
                        cPol.consVar[j,i]+=s.income.mProb[k-2]*(s.income.yDefGrid[i]+s.income.mBounds[1]-cPol.consGrid[j,i])^2
                        cDivYPol.consDivYVar[j,i]+=s.income.mProb[k-2]*((s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[k-2])-cDivYPol.consDivYGrid[j,i])^2
                    end
                end
                #Set the index of the first relevant entry of the borowing policy function
                aPolicyInd=s.pol.firstRepayInd[j,i]

                #Add the contribution of default in the interval in which the threshold lies to the output value
                cPol.consVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*(s.income.yDefGrid[i]+s.income.mBounds[1]-cPol.consGrid[j,i])^2
                cDivYPol.consDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes*((s.income.yDefGrid[i]+s.income.mBounds[1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consDivYGrid[j,i])^2

                tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.defThreshold[j,i]-lastMVal)/s.income.mRes

                #Update the last value of m observed to the threshold level of m at which default occurs
                lastMVal=s.pol.defThreshold[j,i]

                while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
                    if (s.pol.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1]-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*((s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        cPol.consVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1]-cPol.consGrid[j,i])^2
                        cDivYPol.consDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*((s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consDivYGrid[j,i])^2


                        tempPSum+=s.income.mProb[mUBInd-1]*(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                        lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                        aPolicyInd+=1
                    else
                        cPol.consRVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1]-cPol.consRGrid[j,i])^2
                        cDivYPol.consRDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*((s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consRDivYGrid[j,i])^2

                        cPol.consVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1]-cPol.consGrid[j,i])^2
                        cDivYPol.consDivYVar[j,i]+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*((s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]+s.income.mMidPoints[mUBInd-1])/(s.income.yGrid[i]+s.income.mMidPoints[mUBInd-1])-cDivYPol.consDivYGrid[j,i])^2

                        tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                        lastMVal=s.income.mGrid[mUBInd]
                        mUBInd+=1
                    end
                end
                tempPSumInv=tempPSum^(-1)
                cPol.consVar[j,i]*=tempPSumInv
                cDivYPol.consDivYVar[j,i]*=tempPSumInv

                cPol.consRVar[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)
                cDivYPol.consRDivYVar[j,i]*=(1.0-apDefPol.defaultProb[j,i])^(-1)

            else
                cPol.consVar[j,i]=0.0
                cDivYPol.consDivYVar[j,i]=cDivYPol.consDivYDInitVar[i]
            end

        end
    end




    apDefPol.repayProb.=1.0.-apDefPol.defaultProb


    return MeanPolicies(cPol,cDivYPol,apDefPol,uPol,yPol)
end




function writePolicies(s::LongTermBondRDBorrowEval{F,S,T},mPol::MeanPolicies{F,S},filePrefix::String,fileDir::String) where{F<:Real,S<:Integer,T}
    writecsv(
        joinpath(fileDir, filePrefix*"_apRGrid.csv"),
        mPol.apDefPol.apRGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_apRDivYGrid.csv"),
        mPol.apDefPol.apRDivYGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_consRGrid.csv"),
        mPol.cPol.consRGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_consRDivYGrid.csv"),
        mPol.cDivYPol.consRDivYGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_consGrid.csv"),
        mPol.cPol.consGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_consDivYGrid.csv"),
        mPol.cDivYPol.consDivYGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_consDivYDInitGrid.csv"),
        mPol.cDivYPol.consDivYDInitGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_consDivYDFutGrid.csv"),
        mPol.cDivYPol.consDivYDFutGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_uGrid.csv"),
        mPol.uPol.uGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_uDInitGrid.csv"),
        mPol.uPol.uDInitGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_uDFutGrid.csv"),
        mPol.uPol.uDFutGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_defaultProb.csv"),
        mPol.apDefPol.defaultProb
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_repayProb.csv"),
        mPol.apDefPol.repayProb
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_lastAlwaysDefInd.csv"),
        mPol.apDefPol.lastAlwaysDefInd
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_yGrid.csv"),
        s.income.yGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_yDefGrid.csv"),
        s.income.yDefGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_aGrid.csv"),
        s.aGrid
    )
end
