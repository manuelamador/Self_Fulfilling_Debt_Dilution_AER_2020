function genStIncDist(
    m::LongTermBondRDBorrowSpec{F,S},
    s::LongTermBondRDBorrowEval{F,S,T},
    tol::F,
    maxIter::S
) where{F<:Real,S<:Integer,T}

    oldDist=ones(m.yParams.yPoints)/m.yParams.yPoints
    newDist=zeros(m.yParams.yPoints)

    iCount=1
    dDist=1.0+tol

    while (iCount<=maxIter)&&(dDist>tol)
        mul!(newDist,s.income.yTMat,oldDist)
        dDist=maximum(abs.(newDist-oldDist))
        if div(iCount,max(div(maxIter,10),1))==(iCount/max(div(maxIter,10),1))
            # println([iCount,dDist])
            print(".")
        end
        copyto!(oldDist,newDist)
        iCount+=1
    end
    deflFac=sum(newDist)^(-1)
    return newDist.*deflFac
end


function simulatePathsMC(m::LongTermBondRDBorrowSpec{F,S},s::LongTermBondRDBorrowEval{F,S,T},bigT::S,bigN::S,trim::S,trimDef::S) where{F<:Real,S<:Integer,T}
    cSim=zeros(bigT,bigN)
    spreadSim=zeros(bigT,bigN)
    tbSim=zeros(bigT,bigN)
    mvSim=zeros(bigT,bigN)
    yIndSim=zeros(S,bigT,bigN)
    ySim=zeros(bigT,bigN)
    mSim=zeros(bigT,bigN)
    ymSim=zeros(bigT,bigN)
    aIndSim=zeros(S,bigT,bigN)
    aSim=zeros(bigT,bigN)
    apSim=zeros(bigT,bigN)
    apIndSim=zeros(S,bigT,bigN)

    defSim=zeros(Bool,bigT,bigN)
    inDefSim=zeros(Bool,bigT,bigN)
    noDefDuration=zeros(S,bigT,bigN)
    inSample=zeros(Bool,bigT,bigN)
    inDefSample=zeros(Bool,bigT,bigN)

    stDist=genStIncDist(m,s,1e-10,1000)

    stDistCDF=cumsum(stDist)

    mDistTN=Distributions.truncated(
        Normal(
            m.mParams.mu,
            sqrt(m.mParams.epsilon2)
        ),
        s.income.mBounds[1],
        s.income.mBounds[2]
    )
    yTDistCDF=cumsum(s.income.yTMat,dims=1)


    for littleN in 1:bigN
        mLittleT=1
        while mLittleT<=bigT

            mSim[mLittleT,littleN]=rand(mDistTN)
            mLittleT+=1
        end
        for littleT in 1:bigT
            if littleT==1
                yrand=rand()
                yIndSim[1,littleN]=searchsortedfirst(stDistCDF,yrand)
                aIndSim[1,littleN]=s.a0Ind
                aSim[1,littleN]=0.0
            else
                yrand=rand()
                yIndSim[littleT,littleN]=searchsortedfirst(yTDistCDF[:,yIndSim[littleT-1,littleN]],yrand)
                if inDefSim[littleT-1,littleN]==true
                    reenterRand=rand()
                    if reenterRand<m.theta
                        inDefSim[littleT,littleN]=false
                        noDefDuration[littleT,littleN]=1
                        aIndSim[littleT,littleN]=s.a0Ind
                    else
                        inDefSim[littleT,littleN]=true
                        noDefDuration[littleT,littleN]=0
                    end
                end
            end

            if inDefSim[littleT,littleN]==true
                ySim[littleT,littleN]=s.income.yDefGrid[yIndSim[littleT,littleN]]
                if yIndSim[littleT,littleN]!=m.yParams.yPoints
                    ymSim[littleT,littleN]=ySim[littleT,littleN]+mSim[littleT,littleN]
                else

                    ymSim[littleT,littleN]=ySim[littleT,littleN]
                end
                cSim[littleT,littleN]=ymSim[littleT,littleN]

            elseif yIndSim[littleT,littleN]==m.yParams.yPoints
                if s.pol.alwaysDefault[aIndSim[littleT,littleN],yIndSim[littleT,littleN]]==true
                    defSim[littleT,littleN]=true
                    inDefSim[littleT,littleN]=true
                    noDefDuration[littleT,littleN]=0
                    ySim[littleT,littleN]=s.income.yDefGrid[yIndSim[littleT,littleN]]
                    ymSim[littleT,littleN]=ySim[littleT,littleN]
                    cSim[littleT,littleN]=ymSim[littleT,littleN]
                else
                    defSim[littleT,littleN]=false
                    inDefSim[littleT,littleN]=false
                    if littleT!=1
                        noDefDuration[littleT,littleN]=noDefDuration[littleT-1,littleN]+1
                    end
                    ySim[littleT,littleN]=s.income.yGrid[yIndSim[littleT,littleN]]
                    ymSim[littleT,littleN]=ySim[littleT,littleN]

                    apIndSim[littleT,littleN]=s.pol.apHighM[aIndSim[littleT,littleN],yIndSim[littleT,littleN]]
                    apSim[littleT,littleN]=s.aGrid[apIndSim[littleT,littleN]]
                    cSim[littleT,littleN]=s.consM0Grid[apIndSim[littleT,littleN],aIndSim[littleT,littleN],yIndSim[littleT,littleN]]
                    tbSim[littleT,littleN]=ymSim[littleT,littleN]-cSim[littleT,littleN]
                    spreadSim[littleT,littleN]=(m.lambda+(1.0-m.lambda)*m.coup)/s.qGrid[apIndSim[littleT,littleN],yIndSim[littleT,littleN]]-m.lambda
                    mvSim[littleT,littleN]=s.qGrid[apIndSim[littleT,littleN],yIndSim[littleT,littleN]]*apSim[littleT,littleN]
                    if littleT<bigT
                        aIndSim[littleT+1,littleN]=apIndSim[littleT,littleN]
                        aSim[littleT+1,littleN]=apSim[littleT,littleN]
                    end
                end

            elseif (s.pol.alwaysDefault[aIndSim[littleT,littleN],yIndSim[littleT,littleN]]==true)||((s.pol.neverDefault[aIndSim[littleT,littleN],yIndSim[littleT,littleN]]==false)&&(mSim[littleT,littleN]<s.pol.defThreshold[aIndSim[littleT,littleN],yIndSim[littleT,littleN]]))
                defSim[littleT,littleN]=true
                inDefSim[littleT,littleN]=true
                noDefDuration[littleT,littleN]=0
                ySim[littleT,littleN]=s.income.yDefGrid[yIndSim[littleT,littleN]]
                ymSim[littleT,littleN]=ySim[littleT,littleN]+s.income.mBounds[1]
                cSim[littleT,littleN]=ymSim[littleT,littleN]
            else
                defSim[littleT,littleN]=false
                inDefSim[littleT,littleN]=false
                if littleT!=1
                    noDefDuration[littleT,littleN]=noDefDuration[littleT-1,littleN]+1
                end
                ySim[littleT,littleN]=s.income.yGrid[yIndSim[littleT,littleN]]
                ymSim[littleT,littleN]=ySim[littleT,littleN]+mSim[littleT,littleN]
                apPolLength=s.pol.mListLength[aIndSim[littleT,littleN],yIndSim[littleT,littleN]]
                aPolicyInd=searchsortedfirst(s.pol.mAPThreshold[yIndSim[littleT,littleN]][1:apPolLength,aIndSim[littleT,littleN]],mSim[littleT,littleN])
                apIndSim[littleT,littleN]=s.pol.apPolicy[yIndSim[littleT,littleN]][aPolicyInd,aIndSim[littleT,littleN]]
                apSim[littleT,littleN]=s.aGrid[apIndSim[littleT,littleN]]
                cSim[littleT,littleN]=s.consM0Grid[apIndSim[littleT,littleN],aIndSim[littleT,littleN],yIndSim[littleT,littleN]]+mSim[littleT,littleN]
                tbSim[littleT,littleN]=ymSim[littleT,littleN]-cSim[littleT,littleN]
                spreadSim[littleT,littleN]=(m.lambda+(1.0-m.lambda)*m.coup)/s.qGrid[apIndSim[littleT,littleN],yIndSim[littleT,littleN]]-m.lambda
                mvSim[littleT,littleN]=s.qGrid[apIndSim[littleT,littleN],yIndSim[littleT,littleN]]*apSim[littleT,littleN]
                if littleT<bigT
                    aIndSim[littleT+1,littleN]=apIndSim[littleT,littleN]
                    aSim[littleT+1,littleN]=apSim[littleT,littleN]
                end
            end
        end
        for littleT in (trim+1):bigT
            if (noDefDuration[littleT,littleN]>trimDef)
                inSample[littleT,littleN]=true
            end
            if (noDefDuration[littleT-1,littleN]>(trimDef-1))
                inDefSample[littleT,littleN]=true
            end
        end
    end

    meanAPFVDivY=mean(apSim[inSample]./ymSim[inSample])
    meanAPMVDivY=mean(mvSim[inSample]./ymSim[inSample])
    meanAFVNDDivY=mean(aSim[inSample]./ymSim[inSample])
    meanAFVDivY=mean(aSim[inDefSample]./ymSim[inDefSample])

    defRate=mean(defSim[inDefSample])
    meanSpread=mean((1.0 .+spreadSim[inSample]).^4 .-m.R^4)
    volSpread=sqrt(var((1.0 .+spreadSim[inSample]).^4 .-m.R^4))
    volCDivVolY=sqrt(var(log.(cSim[inSample]))/var(log.(ymSim[inSample])))
    volTB=sqrt(var(tbSim[inSample]./ymSim[inSample]))
    corTBLogY=cor(tbSim[inSample]./ymSim[inSample],log.(ymSim[inSample]))

    corSpreadLogY=cor((1.0 .+spreadSim[inSample]).^4 .-m.R^4,log.(ymSim[inSample]))
    corSpreadAPDivYM=cor((1.0 .+spreadSim[inSample]).^4 .-m.R^4,apSim[inSample]./ymSim[inSample])
    corSpreadTB=cor((1.0 .+spreadSim[inSample]).^4 .-m.R^4,tbSim[inSample]./ymSim[inSample])

    inNoRDSample=inSample.*(yIndSim.!=m.yParams.yPoints)
    inNoRDDefSample=inDefSample.*(yIndSim.!=m.yParams.yPoints)

    meanAPFVDivYNRD=mean(apSim[inNoRDSample]./ymSim[inNoRDSample])
    meanAPMVDivYNRD=mean(mvSim[inNoRDSample]./ymSim[inNoRDSample])
    meanAFVNDDivYNRD=mean(aSim[inNoRDSample]./ymSim[inNoRDSample])
    meanAFVDivYNRD=mean(aSim[inNoRDDefSample]./ymSim[inNoRDDefSample])

    defRateNRD=mean(defSim[inNoRDDefSample])
    meanSpreadNRD=mean((1.0 .+spreadSim[inNoRDSample]).^4 .-m.R^4)
    volSpreadNRD=sqrt(var((1.0 .+spreadSim[inNoRDSample]).^4 .-m.R^4))
    volCDivVolYNRD=sqrt(var(log.(cSim[inNoRDSample]))/var(log.(ymSim[inNoRDSample])))
    volTBNRD=sqrt(var(tbSim[inNoRDSample]./ymSim[inNoRDSample]))
    corTBLogYNRD=cor(tbSim[inNoRDSample]./ymSim[inNoRDSample],log.(ymSim[inNoRDSample]))

    corSpreadLogYNRD=cor((1.0 .+spreadSim[inNoRDSample]).^4 .-m.R^4,log.(ymSim[inNoRDSample]))
    corSpreadAPDivYMNRD=cor((1.0 .+spreadSim[inNoRDSample]).^4 .-m.R^4,apSim[inNoRDSample]./ymSim[inNoRDSample])
    corSpreadTBNRD=cor((1.0 .+spreadSim[inNoRDSample]).^4 .-m.R^4,tbSim[inNoRDSample]./ymSim[inNoRDSample])

    momentValsFull=[meanAPFVDivY,meanAPMVDivY,meanAFVNDDivY,meanAFVDivY,defRate,meanSpread,volSpread,volCDivVolY,volTB,corTBLogY,corSpreadLogY,corSpreadAPDivYM,corSpreadTB]
    momentValsFullNRD=[meanAPFVDivYNRD,meanAPMVDivYNRD,meanAFVNDDivYNRD,meanAFVDivYNRD,defRateNRD,meanSpreadNRD,volSpreadNRD,volCDivVolYNRD,volTBNRD,corTBLogYNRD,corSpreadLogYNRD,corSpreadAPDivYMNRD,corSpreadTBNRD]

    MCPathsOut=MCPaths(cSim,ySim,mSim,ymSim,aSim,apSim,tbSim,spreadSim,mvSim,defSim,inDefSim,yIndSim,aIndSim,apIndSim,noDefDuration,inSample,inDefSample)
    return MCPathsOut,momentValsFull,momentValsFullNRD
end


function simulatePathsMC(m::LongTermBondRDBorrowSpec{F,S},s::LongTermBondRDBorrowEval{F,S,T},bigT::S,bigN::S,trim::S,trimDef::S,thetaAlt::F) where{F<:Real,S<:Integer,T}
    mNew=LongTermBondRDBorrowSpec(m.beta,thetaAlt,m.gamma,m.R,m.lambda,m.coup,m.hpen0,m.hpen1,m.aBounds,m.aPoints,m.yParams,m.mParams,m.simplePen,m.mixFacQ,m.mixFacV)

    return simulatePathsMC(mNew,s,bigT,bigN,trim,trimDef)
end
