function makeVDFuncSave(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},yDefGrid::Array{F,1},theta::F) where{F<:Real,S<:Integer,T}
    ydim=m.yParams.yPoints
    adim=m.aPoints

    vDInitFlow=zeros(ydim)
    vDFutFlow=zeros(ydim)
    vDInitGrid=zeros(ydim)
    vDFutGrid=zeros(ydim)
    EVDGrid=zeros(ydim)

    for i in 1:(ydim-1)
        vDFutFlow[i]=0.0
        for k in 1:(m.mParams.mPoints-1)
            vDFutFlow[i]+=s.income.mProb[k]*(log(yDefGrid[i]+s.income.mGrid[k])-log(yDefGrid[i]+s.income.mGrid[k+1]))/s.income.mRes
        end
        vDInitFlow[i]=u(m,yDefGrid[i]+s.income.mBounds[1])

    end

    vDFutFlow[ydim]=u(m,yDefGrid[ydim])
    vDInitFlow[ydim]=u(m,yDefGrid[ydim])

    vDFutSolnMat=inv(Matrix{F}(I,ydim,ydim).-m.beta.*(1.0-theta).*s.income.yTMat')

    vDFutGrid.=vDFutSolnMat*(vDFutFlow.+(m.beta*theta).*s.VF.EVGrid[s.a0Ind,:])
    EVDGrid.=(1.0-theta).*s.income.yTMat'*vDFutGrid.+theta.*s.VF.EVGrid[s.a0Ind,:]
    vDInitGrid.=vDInitFlow.+m.beta*EVDGrid

    tempVLow=zeros(ydim)
    for i in 1:ydim-1
        tempVLow[i]=u(m,s.consM0Grid[s.pol.apLowM[1,i],1,i]+s.income.mBounds[1])+m.beta*s.VF.EVGrid[s.pol.apLowM[1,i],i]
    end
    tempVLow[ydim]=u(m,s.consM0Grid[s.pol.apLowM[1,ydim],1,ydim])+m.beta*s.VF.EVGrid[s.pol.apLowM[1,ydim],ydim]

    print(".")
    # println(maximum(abs.(vDFutGrid.-(vDFutFlow.+m.beta*EVDGrid))))
    # println(maximum(vDFutGrid .- s.VF.vGrid[1,:]))
    # println(maximum(vDInitGrid .- s.VF.vGrid[1,:]))
    # println(maximum(vDFutGrid[1:ydim-1] .- s.VF.vGrid[1,1:ydim-1]))
    # println(maximum(vDInitGrid[1:ydim-1] .- s.VF.vGrid[1,1:ydim-1]))
    # println(maximum(vDInitGrid[1:ydim-1] .- tempVLow[1:ydim-1]))

    return VDFuncSave(vDInitFlow,vDFutFlow,vDInitGrid,vDFutGrid,EVDGrid)
end


function setupEqmCompletion(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},aExtBounds::Array{F,1},aExtPoints::S,theta::F,hpen0::F,hpen1::F) where{F<:Real,S<:Integer,T}
    ydim=m.yParams.yPoints
    adim=m.aPoints

    @assert aExtBounds[2]==m.aBounds[1]
    @assert aExtBounds[1]<aExtBounds[2]

    yDefGrid=(1.0 .- max.(0.0,hpen0 .+ hpen1 .* s.income.yGrid)).*s.income.yGrid

    @assert minimum(yDefGrid)>0.0
    @assert minimum(diff(yDefGrid[1:end-1]))>0.0

    aExtRes=(aExtBounds[2]-aExtBounds[1])/aExtPoints
    aGridExtension=collect(LinRange(aExtBounds[1],aExtBounds[2]-aExtRes,aExtPoints))

    qGrid=ones(F,adim+aExtPoints,ydim)*s.qBase
    qGrid[1:aExtPoints,:].=0.0
    aGrid=vcat(aGridExtension,s.aGrid)
    aGridIncr=aGrid.-(1.0-m.lambda)*aGrid'

    netRevM0A0Grid=zeros(F,adim+aExtPoints,ydim)
    consM0Grid=zeros(F,adim+aExtPoints,adim+aExtPoints,ydim)

    for i in 1:ydim
        for j in 1:(adim+aExtPoints)
            netRevM0A0Grid[j,i]=s.income.yGrid[i].+(m.lambda+(1.0-m.lambda)*m.coup)*aGrid[j]
            view(consM0Grid,:,j,i).=netRevM0A0Grid[j,i].-view(qGrid,:,i).*view(aGridIncr,:,j)

        end
    end

    vExtension=zeros(aExtPoints,ydim)
    EVExtension=zeros(aExtPoints,ydim)

    for j in 1:aExtPoints
        vExtension[j,:].=s.VF.vGrid[1,:].+(s.VF.vGrid[2,:].-s.VF.vGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*10.0
        EVExtension[j,:].=s.VF.EVGrid[1,:].+(s.VF.EVGrid[2,:].-s.VF.EVGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*10.0
    end
    vExtension[1:aExtPoints,m.yParams.yPoints].=s.VF.vGrid[1,m.yParams.yPoints]

    VFExtension=VFuncSave(vExtension,EVExtension)
    vFull=vcat(vExtension,s.VF.vGrid)
    EVFull=vcat(EVExtension,s.VF.EVGrid)
    VFFull=VFuncSave(vFull,EVFull)

    VFD=makeVDFuncSave(m,s,yDefGrid,theta)

    apLowMExtension=ones(S,aExtPoints,ydim)
    apHighMExtension=ones(S,aExtPoints,ydim)
    mListLengthExtension=zeros(S,aExtPoints,ydim)
    mAPThresholdExtension=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyExtension=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFExtension=sparse(zeros(F,adim+aExtPoints,aExtPoints))
    sDummySExtension=sparse(zeros(S,adim+aExtPoints,aExtPoints))
    for i in 1:ydim
        mAPThresholdExtension[i]=deepcopy(sDummyFExtension)
        apPolicyExtension[i]=deepcopy(sDummySExtension)
    end

    qSumExtension=deepcopy(qGrid[1:aExtPoints,:])

    polExtension=PoliciesSave(apLowMExtension,apHighMExtension,mListLengthExtension,mAPThresholdExtension,apPolicyExtension)

    apLowMFull=ones(S,adim+aExtPoints,ydim)
    apHighMFull=ones(S,adim+aExtPoints,ydim)
    mListLengthFull=zeros(S,adim+aExtPoints,ydim)
    mAPThresholdFull=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyFull=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFFull=sparse(zeros(F,adim+aExtPoints,adim+aExtPoints))
    sDummySFull=sparse(zeros(S,adim+aExtPoints,adim+aExtPoints))
    for i in 1:ydim
        mAPThresholdFull[i]=deepcopy(sDummyFFull)
        apPolicyFull[i]=deepcopy(sDummySFull)
    end

    apLowMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apLowM.+aExtPoints
    apHighMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apHighM.+aExtPoints
    mListLengthFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.mListLength

    for i in 1:ydim
        for j in (aExtPoints+1):(adim+aExtPoints)
            for jj in 1:s.pol.mListLength[j-aExtPoints,i]
                mAPThresholdFull[i][jj,j]=s.pol.mAPThreshold[i][jj,j-aExtPoints]
                apPolicyFull[i][jj,j]=s.pol.apPolicy[i][jj,j-aExtPoints]+aExtPoints
            end
        end
    end

    polFull=PoliciesSave(apLowMFull,apHighMFull,mListLengthFull,mAPThresholdFull,apPolicyFull)

    defThreshold=ones(F,adim+aExtPoints,ydim)*(m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    firstRepayInd=ones(S,adim+aExtPoints,ydim)
    alwaysDefault=zeros(Bool,adim+aExtPoints,ydim)
    neverDefault=ones(Bool,adim+aExtPoints,ydim)

    defThreshold[1:aExtPoints,ydim].=(-m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    alwaysDefault[1:aExtPoints,ydim].=true
    neverDefault[1:aExtPoints,ydim].=false

    attemptedSqRootNeg=zeros(Bool,aExtPoints,ydim)

    dPol=DPoliciesSave(defThreshold,firstRepayInd,alwaysDefault,neverDefault,attemptedSqRootNeg)

    apSearchParams=genBisectSearchGrids(aExtPoints)

    return LongTermBondRDSaveExtension(yDefGrid,aExtPoints,aGrid,aGridIncr,qGrid,netRevM0A0Grid,consM0Grid,VFExtension,VFFull,VFD,polExtension,qSumExtension,polFull,dPol,apSearchParams,[true])
end




function setupEqmCompletion(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},aGridExtension::Array{F,1},theta::F,hpen0::F,hpen1::F) where{F<:Real,S<:Integer,T}
    ydim=m.yParams.yPoints
    adim=m.aPoints

    @assert aGridExtension[end]==m.aBounds[1]
    @assert issorted(aGridExtension)

    yDefGrid=(1.0 .- max.(0.0,hpen0 .+ hpen1 .* s.income.yGrid)).*s.income.yGrid

    aExtPoints=length(aGridExtension)-1


    @assert minimum(yDefGrid)>0.0
    @assert minimum(diff(yDefGrid[1:end-1]))>0.0

    qGrid=ones(F,adim+aExtPoints,ydim)*s.qBase
    qGrid[1:aExtPoints,:]=0.0
    aGrid=vcat(aGridExtension[1:end-1],s.aGrid)
    aGridIncr=aGrid.-(1.0-m.lambda)*aGrid'

    netRevM0A0Grid=zeros(F,adim+aExtPoints,ydim)
    consM0Grid=zeros(F,adim+aExtPoints,adim+aExtPoints,ydim)

    for i in 1:ydim
        for j in 1:(adim+aExtPoints)
            netRevM0A0Grid[j,i]=s.income.yGrid[i].+(m.lambda+(1.0-m.lambda)*m.coup)*aGrid[j]
            view(consM0Grid,:,j,i).=netRevM0A0Grid[j,i].-view(qGrid,:,i).*view(aGridIncr,:,j)

        end
    end

    VFD=makeVDFuncSave(m,s,yDefGrid,theta)

    vExtension=zeros(aExtPoints,ydim)
    EVExtension=zeros(aExtPoints,ydim)

    for j in 1:aExtPoints
        vExtension[j,:].=s.VF.vGrid[1,:].+(s.VF.vGrid[2,:].-s.VF.vGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*10.0
        EVExtension[j,:].=s.VF.EVGrid[1,:].+(s.VF.EVGrid[2,:].-s.VF.EVGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*10.0
    end

    vExtension[1:aExtPoints,m.yParams.yPoints].=VFD.vDInitGrid[m.yParams.yPoints]

    VFExtension=VFuncSave(vExtension,EVExtension)
    vFull=vcat(vExtension,s.VF.vGrid)
    EVFull=vcat(EVExtension,s.VF.EVGrid)
    VFFull=VFuncSave(vFull,EVFull)

    apLowMExtension=ones(S,aExtPoints,ydim)
    apHighMExtension=ones(S,aExtPoints,ydim)
    mListLengthExtension=zeros(S,aExtPoints,ydim)
    mAPThresholdExtension=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyExtension=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFExtension=sparse(zeros(F,adim+aExtPoints,aExtPoints))
    sDummySExtension=sparse(zeros(S,adim+aExtPoints,aExtPoints))
    for i in 1:ydim
        mAPThresholdExtension[i]=deepcopy(sDummyFExtension)
        apPolicyExtension[i]=deepcopy(sDummySExtension)
    end
    qSumExtension=deepcopy(qGrid[1:aExtPoints,:])

    polExtension=PoliciesSave(apLowMExtension,apHighMExtension,mListLengthExtension,mAPThresholdExtension,apPolicyExtension)

    apLowMFull=ones(S,adim+aExtPoints,ydim)
    apHighMFull=ones(S,adim+aExtPoints,ydim)
    mListLengthFull=zeros(S,adim+aExtPoints,ydim)
    mAPThresholdFull=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyFull=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFFull=sparse(zeros(F,adim+aExtPoints,adim+aExtPoints))
    sDummySFull=sparse(zeros(S,adim+aExtPoints,adim+aExtPoints))
    for i in 1:ydim
        mAPThresholdFull[i]=deepcopy(sDummyFFull)
        apPolicyFull[i]=deepcopy(sDummySFull)
    end

    apLowMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apLowM.+aExtPoints
    apHighMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apHighM.+aExtPoints
    mListLengthFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.mListLength

    for i in 1:ydim
        for j in (aExtPoints+1):(adim+aExtPoints)
            for jj in 1:s.pol.mListLength[j-aExtPoints,i]
                mAPThresholdFull[i][jj,j]=s.pol.mAPThreshold[i][jj,j-aExtPoints]
                apPolicyFull[i][jj,j]=s.pol.apPolicy[i][jj,j-aExtPoints]+aExtPoints
            end
        end
    end

    polFull=PoliciesSave(apLowMFull,apHighMFull,mListLengthFull,mAPThresholdFull,apPolicyFull)

    defThreshold=ones(F,adim+aExtPoints,ydim)*(m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    firstRepayInd=ones(S,adim+aExtPoints,ydim)
    alwaysDefault=zeros(Bool,adim+aExtPoints,ydim)
    neverDefault=ones(Bool,adim+aExtPoints,ydim)

    defThreshold[1:aExtPoints,ydim].=(-m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    alwaysDefault[1:aExtPoints,ydim].=true
    neverDefault[1:aExtPoints,ydim].=false

    attemptedSqRootNeg=zeros(Bool,aExtPoints,ydim)

    dPol=DPoliciesSave(defThreshold,firstRepayInd,alwaysDefault,neverDefault,attemptedSqRootNeg)
    apSearchParams=genBisectSearchGrids(aExtPoints)

    return LongTermBondRDSaveExtension(yDefGrid,aExtPoints,aGrid,aGridIncr,qGrid,netRevM0A0Grid,consM0Grid,VFExtension,VFFull,VFD,polExtension,qSumExtension,polFull,dPol,apSearchParams,[true])
end



function solveYDefRDSave(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},theta::F,hpen0::F,hpen1::F,shiftRDVal::F) where{F<:Real,S<:Integer,T}
    ydim=m.yParams.yPoints
    adim=m.aPoints

    vDInitFlow=zeros(ydim)
    vDFutFlow=zeros(ydim)
    vDInitGrid=zeros(ydim)
    vDFutGrid=zeros(ydim)
    EVDGrid=zeros(ydim)

    yDefGrid=zeros(ydim)

    for i in 1:(ydim-1)
        yDefGrid[i]=s.income.yGrid[i]-max(0,hpen0*s.income.yGrid[i]+hpen1*s.income.yGrid[i]^2)
        vDFutFlow[i]=0.0
        for k in 1:(m.mParams.mPoints-1)
            vDFutFlow[i]+=s.income.mProb[k]*(log(yDefGrid[i]+s.income.mGrid[k])-log(yDefGrid[i]+s.income.mGrid[k+1]))/s.income.mRes
        end
        vDInitFlow[i]=u(m,yDefGrid[i]+s.income.mBounds[1])
    end
    vDFutNRDSolnMat=inv(Matrix{F}(I,ydim-1,ydim-1).-m.beta.*(1.0-theta).*s.income.yTMat[1:ydim-1,1:ydim-1]')
    vDFutFlowHat=vDFutFlow[1:ydim-1].+s.income.yTMat[ydim,1:ydim-1].*m.beta*(1.0-theta)*(s.VF.vGrid[1,ydim]-shiftRDVal)
    vDFutGridHat=vDFutNRDSolnMat*(vDFutFlowHat.+(m.beta*theta).*s.VF.EVGrid[s.a0Ind,1:ydim-1])
    vDFutGrid[1:ydim-1].=vDFutGridHat
    vDFutGrid[ydim]=(s.VF.vGrid[1,ydim]-shiftRDVal)
    EVDGrid.=(1.0-theta).*s.income.yTMat'*vDFutGrid.+theta.*s.VF.EVGrid[s.a0Ind,:]

    vDFutFlow[ydim]=vDFutGrid[ydim]-m.beta*EVDGrid[ydim]
    yDefGrid[ydim]=uInverse(m,vDFutFlow[ydim])
    vDInitFlow[ydim]=vDFutFlow[ydim]

    vDInitGrid.=vDInitFlow.+m.beta*EVDGrid

    tempVLow=zeros(ydim)
    for i in 1:ydim-1
        tempVLow[i]=u(m,s.consM0Grid[s.pol.apLowM[1,i],1,i]+s.income.mBounds[1])+m.beta*s.VF.EVGrid[s.pol.apLowM[1,i],i]
    end
    tempVLow[ydim]=u(m,s.consM0Grid[s.pol.apLowM[1,ydim],1,ydim])+m.beta*s.VF.EVGrid[s.pol.apLowM[1,ydim],ydim]

    # println("\n", minimum(tempVLow[1:end-1].-vDInitGrid[1:end-1]))
    # println(minimum(tempVLow.-vDInitGrid))

    return VDFuncSave(vDInitFlow,vDFutFlow,vDInitGrid,vDFutGrid,EVDGrid),yDefGrid
end


function setupEqmCompletionYDefRDImpl(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},aExtBounds::Array{F,1},aExtPoints::S,theta::F,hpen0::F,hpen1::F) where{F<:Real,S<:Integer,T}
    ydim=m.yParams.yPoints
    adim=m.aPoints

    @assert aExtBounds[2]==m.aBounds[1]
    @assert aExtBounds[1]<aExtBounds[2]

    aExtRes=(aExtBounds[2]-aExtBounds[1])/aExtPoints
    aGridExtension=collect(LinRange(aExtBounds[1],aExtBounds[2]-aExtRes,aExtPoints))
    aGrid=vcat(aGridExtension,s.aGrid)
    aGridIncr=aGrid.-(1.0-m.lambda)*aGrid'

    qGrid=ones(F,adim+aExtPoints,ydim)*s.qBase
    qGrid[1:aExtPoints,:].=0.0

    netRevM0A0Grid=zeros(F,adim+aExtPoints,ydim)
    consM0Grid=zeros(F,adim+aExtPoints,adim+aExtPoints,ydim)

    for i in 1:ydim
        for j in 1:(adim+aExtPoints)
            netRevM0A0Grid[j,i]=s.income.yGrid[i].+(m.lambda+(1.0-m.lambda)*m.coup)*aGrid[j]
            view(consM0Grid,:,j,i).=netRevM0A0Grid[j,i].-view(qGrid,:,i).*view(aGridIncr,:,j)
        end
    end

    vExtension=zeros(aExtPoints,ydim)
    EVExtension=zeros(aExtPoints,ydim)

    for j in 1:aExtPoints
        vExtension[j,:].=s.VF.vGrid[1,:].+(s.VF.vGrid[2,:].-s.VF.vGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*10.0
        EVExtension[j,:].=s.VF.EVGrid[1,:].+(s.VF.EVGrid[2,:].-s.VF.EVGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*10.0
    end

    VFD,yDefGrid=solveYDefRDSave(m,s,theta,hpen0,hpen1,1e-13)

    vExtension[1:aExtPoints,m.yParams.yPoints].=VFD.vDInitGrid[m.yParams.yPoints]

    VFExtension=VFuncSave(vExtension,EVExtension)
    vFull=vcat(vExtension,s.VF.vGrid)
    EVFull=vcat(EVExtension,s.VF.EVGrid)
    VFFull=VFuncSave(vFull,EVFull)

    apLowMExtension=ones(S,aExtPoints,ydim)
    apHighMExtension=ones(S,aExtPoints,ydim)
    mListLengthExtension=zeros(S,aExtPoints,ydim)
    mAPThresholdExtension=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyExtension=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFExtension=sparse(zeros(F,adim+aExtPoints,aExtPoints))
    sDummySExtension=sparse(zeros(S,adim+aExtPoints,aExtPoints))
    for i in 1:ydim
        mAPThresholdExtension[i]=deepcopy(sDummyFExtension)
        apPolicyExtension[i]=deepcopy(sDummySExtension)
    end
    qSumExtension=deepcopy(qGrid[1:aExtPoints,:])

    polExtension=PoliciesSave(apLowMExtension,apHighMExtension,mListLengthExtension,mAPThresholdExtension,apPolicyExtension)

    apLowMFull=ones(S,adim+aExtPoints,ydim)
    apHighMFull=ones(S,adim+aExtPoints,ydim)
    mListLengthFull=zeros(S,adim+aExtPoints,ydim)
    mAPThresholdFull=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyFull=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFFull=sparse(zeros(F,adim+aExtPoints,adim+aExtPoints))
    sDummySFull=sparse(zeros(S,adim+aExtPoints,adim+aExtPoints))
    for i in 1:ydim
        mAPThresholdFull[i]=deepcopy(sDummyFFull)
        apPolicyFull[i]=deepcopy(sDummySFull)
    end

    apLowMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apLowM.+aExtPoints
    apHighMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apHighM.+aExtPoints
    mListLengthFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.mListLength

    for i in 1:ydim
        for j in (aExtPoints+1):(adim+aExtPoints)
            for jj in 1:s.pol.mListLength[j-aExtPoints,i]
                mAPThresholdFull[i][jj,j]=s.pol.mAPThreshold[i][jj,j-aExtPoints]
                apPolicyFull[i][jj,j]=s.pol.apPolicy[i][jj,j-aExtPoints]+aExtPoints
            end
        end
    end

    polFull=PoliciesSave(apLowMFull,apHighMFull,mListLengthFull,mAPThresholdFull,apPolicyFull)

    defThreshold=ones(F,adim+aExtPoints,ydim)*(m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    firstRepayInd=ones(S,adim+aExtPoints,ydim)
    alwaysDefault=zeros(Bool,adim+aExtPoints,ydim)
    neverDefault=ones(Bool,adim+aExtPoints,ydim)

    defThreshold[1:aExtPoints,ydim].=(-m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    alwaysDefault[1:aExtPoints,ydim].=true
    neverDefault[1:aExtPoints,ydim].=false

    attemptedSqRootNeg=zeros(Bool,aExtPoints,ydim)

    dPol=DPoliciesSave(defThreshold,firstRepayInd,alwaysDefault,neverDefault,attemptedSqRootNeg)

    apSearchParams=genBisectSearchGrids(aExtPoints)

    return LongTermBondRDSaveExtension(yDefGrid,aExtPoints,aGrid,aGridIncr,qGrid,netRevM0A0Grid,consM0Grid,VFExtension,VFFull,VFD,polExtension,qSumExtension,polFull,dPol,apSearchParams,[true])
end


function setupEqmCompletionYDefRDImpl(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},aGridExtension::Array{F,1},theta::F,hpen0::F,hpen1::F) where{F<:Real,S<:Integer,T}
    ydim=m.yParams.yPoints
    adim=m.aPoints

    @assert aGridExtension[end]==m.aBounds[1]
    @assert issorted(aGridExtension)

    aExtPoints=length(aGridExtension)-1

    aGrid=vcat(aGridExtension[1:end-1],s.aGrid)
    aGridIncr=aGrid.-(1.0-m.lambda)*aGrid'

    qGrid=ones(F,adim+aExtPoints,ydim)*s.qBase
    qGrid[1:aExtPoints,:].=0.0

    netRevM0A0Grid=zeros(F,adim+aExtPoints,ydim)
    consM0Grid=zeros(F,adim+aExtPoints,adim+aExtPoints,ydim)

    for i in 1:ydim
        for j in 1:(adim+aExtPoints)
            netRevM0A0Grid[j,i]=s.income.yGrid[i].+(m.lambda+(1.0-m.lambda)*m.coup)*aGrid[j]
            view(consM0Grid,:,j,i).=netRevM0A0Grid[j,i].-view(qGrid,:,i).*view(aGridIncr,:,j)
        end
    end

    vExtension=zeros(aExtPoints,ydim)
    EVExtension=zeros(aExtPoints,ydim)

    for j in 1:aExtPoints
        vExtension[j,:].=s.VF.vGrid[1,:].+(s.VF.vGrid[2,:].-s.VF.vGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*2.0#10.0
        EVExtension[j,:].=s.VF.EVGrid[1,:].+(s.VF.EVGrid[2,:].-s.VF.EVGrid[1,:])/(s.aGrid[2]-s.aGrid[1])*(aGrid[j]-aGrid[aExtPoints+1])*2.0#10.0
    end
    VFD,yDefGrid=solveYDefRDSave(m,s,theta,hpen0,hpen1,1e-13)

    vExtension[1:aExtPoints,m.yParams.yPoints].=VFD.vDInitGrid[m.yParams.yPoints]

    VFExtension=VFuncSave(vExtension,EVExtension)
    vFull=vcat(vExtension,s.VF.vGrid)
    EVFull=vcat(EVExtension,s.VF.EVGrid)
    VFFull=VFuncSave(vFull,EVFull)


    apLowMExtension=ones(S,aExtPoints,ydim)
    apHighMExtension=ones(S,aExtPoints,ydim)
    mListLengthExtension=zeros(S,aExtPoints,ydim)
    mAPThresholdExtension=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyExtension=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFExtension=sparse(zeros(F,adim+aExtPoints,aExtPoints))
    sDummySExtension=sparse(zeros(S,adim+aExtPoints,aExtPoints))
    for i in 1:ydim
        mAPThresholdExtension[i]=deepcopy(sDummyFExtension)
        apPolicyExtension[i]=deepcopy(sDummySExtension)
    end
    qSumExtension=deepcopy(qGrid[1:aExtPoints,:])

    polExtension=PoliciesSave(apLowMExtension,apHighMExtension,mListLengthExtension,mAPThresholdExtension,apPolicyExtension)

    apLowMFull=ones(S,adim+aExtPoints,ydim)
    apHighMFull=ones(S,adim+aExtPoints,ydim)
    mListLengthFull=zeros(S,adim+aExtPoints,ydim)
    mAPThresholdFull=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicyFull=Array{SparseMatrixCSC{S,S},1}(undef,ydim)

    sDummyFFull=sparse(zeros(F,adim+aExtPoints,adim+aExtPoints))
    sDummySFull=sparse(zeros(S,adim+aExtPoints,adim+aExtPoints))
    for i in 1:ydim
        mAPThresholdFull[i]=deepcopy(sDummyFFull)
        apPolicyFull[i]=deepcopy(sDummySFull)
    end

    apLowMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apLowM.+aExtPoints
    apHighMFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.apHighM.+aExtPoints
    mListLengthFull[(aExtPoints+1):(adim+aExtPoints),:].=s.pol.mListLength

    for i in 1:ydim
        for j in (aExtPoints+1):(adim+aExtPoints)
            for jj in 1:s.pol.mListLength[j-aExtPoints,i]
                mAPThresholdFull[i][jj,j]=s.pol.mAPThreshold[i][jj,j-aExtPoints]
                apPolicyFull[i][jj,j]=s.pol.apPolicy[i][jj,j-aExtPoints]+aExtPoints
            end
        end
    end

    polFull=PoliciesSave(apLowMFull,apHighMFull,mListLengthFull,mAPThresholdFull,apPolicyFull)

    defThreshold=ones(F,adim+aExtPoints,ydim)*(m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    firstRepayInd=ones(S,adim+aExtPoints,ydim)
    alwaysDefault=zeros(Bool,adim+aExtPoints,ydim)
    neverDefault=ones(Bool,adim+aExtPoints,ydim)

    defThreshold[1:aExtPoints,ydim].=(-m.mParams.stdSpan*sqrt(m.mParams.epsilon2)+m.mParams.mu)
    alwaysDefault[1:aExtPoints,ydim].=true
    neverDefault[1:aExtPoints,ydim].=false

    attemptedSqRootNeg=zeros(Bool,aExtPoints,ydim)

    dPol=DPoliciesSave(defThreshold,firstRepayInd,alwaysDefault,neverDefault,attemptedSqRootNeg)

    apSearchParams=genBisectSearchGrids(aExtPoints)

    return LongTermBondRDSaveExtension(yDefGrid,aExtPoints,aGrid,aGridIncr,qGrid,netRevM0A0Grid,consM0Grid,VFExtension,VFFull,VFD,polExtension,qSumExtension,polFull,dPol,apSearchParams,[true])
end


#The function makeLTBUpdate constructs a combinedIncomeProcessSave object based on the contents of a model specification and eval object.
#Its arguments are:
#1. m: a model specification
#2. s: the collection of objects used to solve the model

#Its output is simply a longTermBondRDSaveUpdate object

function makeLTBExtensionUpdate(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S}) where{F<:Real,S<:Integer,T}
    #Set some aliases
    ydim=m.yParams.yPoints
    adim=m.aPoints

    VFExtension=deepcopy(sExt.VFExtension)
    VF=deepcopy(sExt.VF)

    qGridExtension=sExt.qGrid[1:sExt.aExtPoints,:]
    qSumExtension=deepcopy(sExt.qSumExtension)

    maxAlwaysDefInd=ones(S,ydim)
    maxAlwaysDefInd[ydim]=sExt.aExtPoints

    solveMarkH=zeros(Bool,sExt.aExtPoints,ydim)
    solveMarkL=zeros(Bool,sExt.aExtPoints,ydim)
    feasibleSolutionH=zeros(Bool,sExt.aExtPoints,ydim)
    feasibleSolutionL=zeros(Bool,sExt.aExtPoints,ydim)

    return LongTermBondRDExtensionUpdate(VFExtension,VF,qGridExtension,qSumExtension,maxAlwaysDefInd,solveMarkH,solveMarkL,feasibleSolutionH,feasibleSolutionL)
end


function c(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},i::S,j::S,apInd::S,x::U,dMark::Bool) where{F<:Real,U<:Real,S<:Integer,T}
    if dMark==false
        return sExt.netRevM0A0Grid[j,i]-sExt.qGrid[apInd,i]*sExt.aGridIncr[apInd,j]+x
    else
        return sExt.yDefGrid[i]+x
    end
end


function c(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},i::S,j::S,apInd::S,x::U) where{F<:Real,U<:Real,S<:Integer,T}

    return sExt.netRevM0A0Grid[j,i]-sExt.qGrid[apInd,i]*sExt.aGridIncr[apInd,j]+x

end


function c(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},i::S,j::S,apInd::S) where{F<:Real,S<:Integer,T}

    return sExt.netRevM0A0Grid[j,i]-sExt.qGrid[apInd,i]*sExt.aGridIncr[apInd,j]

end


#solveRepayChoice is a function which solve the government's problem under repayment for a specific value of persistent income, transient income, and incoming asset level.
#Its arguments are:
#1. m: the model specification
#2. s: the list of objects used to solve the model
#3. i: the index of the persistent income state
#4. j: the index of the incoming asset state
#5. x: the level of the m shock
#6. lb: the minimum index in the borrowing grid for the grid search
#7. ub: the maximum index in the borowing grid for the grid search

#This function assumes that lb<=ub so the check lb==ub is equivalent to lb<=ub. In the current implementation, it is never passed values which disobey this.
function solveRepayChoice(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},i::S,j::S,x::F,lb::S,ub::S) where{F<:Real,S<:Integer,T}
    #If the lower and upper bound are the same, return the that index and the corresponding value as the solution
    tempC=0.0
    maxC=0.0


    if lb==ub
        tempC=c(m,s,sExt,i,j,lb,x)
        return lb, u(m,tempC)+m.beta*sExt.VF.EVGrid[lb,i],tempC>zero(F)
    else
        #Generate temporary variables containing (in order):
        #1. the maximum consumption observed over levels of next period borrowing examined
        #2. the value of consumption for the level of next period borrowing currently being considered
        #3. the maximum available value thusfar observed for the government
        #4. the index of next period borrowing which results in the value of 3.

        #Initialize these to correspond to the values when the government borrows at the level indicated by lb
        maxC=c(m,s,sExt,i,j,ub,x)
        tempC=c(m,s,sExt,i,j,ub,x)
        maxVR=u(m,tempC)+m.beta*sExt.VF.EVGrid[ub,i]
        maxAPInd=ub

        #Iterate over the remaining values of next period borrowing
        for jj in (ub-1):-1:(lb)

            #Set the value of consumption for the level of next period borrowing currently under consideration
            tempC=c(m,s,sExt,i,j,jj,x)

            #If that value is higher than the maximum value observed thusfar, calculate the government's value under this choice
            if tempC>maxC
                tempVR=u(m,tempC)+m.beta*sExt.VF.EVGrid[jj,i]

                #If the value of the current choice is strictly greater than the highest value thusfar observed, set the maximum value and maximizing value variables accordingly
                if tempVR>maxVR
                    maxVR=tempVR
                    maxAPInd=jj

                end

                #Track that the maximum value of consumption observed thusfar has changed
                maxC=tempC
            end

        end
        tempC=c(m,s,sExt,i,j,maxAPInd,x)

        #Return the maximizing value of the next period borrowing index, the value to the government associated with it,
        #and whether the optimal choice of asset resulted in strictly positive consumption
        return maxAPInd, maxVR,tempC>zero(F)
    end
end



#solveRepayInterior! constructs the full borrowing policy functions and
#default policy functions of the government (should the government repay under any realization of m).
#Its arguments are:
#1. m: the model specification
#2. s: the list of objects used to solve the model
#3. i: the index of the persistent income state
#4. j: the index of the incoming asset state
#5. setExactP: a Boolean variable indicating whether to calculate the exact transition probabilities implied by the policy functions

function solveRepayInterior!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},i::S,j::S) where{F<:Real,S<:Integer,T}

    #Set aliases for the index of assets chosen when m takes its highest possible value and the index of assets chosen when m takes its highest possible value
    lb=sExt.polExtension.apLowM[j,i]
    ub=sExt.polExtension.apHighM[j,i]

    #THIS SHOULD ONLY EVERY BE INVOKED WHEN lb==ub. THE IMPLEMENTATION WHICH PRECEDES THIS IN THE MAIN FUNCTIONS SHOULD BE SUCH THAT
    #lb<=ub WHENEVER solveRepayInterior! IS CALLED.
    #IT ALSO ASSUMES THAT IT IS NEVER CALLED WHEN neverDefault[j,i] IS TRUE. FOR THE CURRENT IMPLEMENTATION WHICH PRECEDES THIS IN THE
    #MAIN FUNCTIONS, THIS IS CERTAINLY ALWAYS THE CASE.

    #If the government chooses the same value at the lowest realization of m that it does at the highest realization of m, use a special method
    #for filling in the contents of the government policy functions.
    if lb>=ub
        #In every case, we set:
        #1. the m threshold at which the government ceases to use the first element of its policy function to the maximum value of m
        #2. the borrowing policy for the first element to be the lower bound
        if sExt.dPol.neverDefault[j,i]==true
            sExt.polExtension.mAPThreshold[i][1,j]=s.income.mBounds[2]
            sExt.polExtension.apPolicy[i][1,j]=lb
            #When the government never defaults in income/debt state (i,j), further set:
            #3. the first index of the borrowing policy function which is relevant for the integration to 1
            #4. the threshold for default to the lowest possible value of the m shock
            #5. the length of the list of thresholds for this state to 1
            sExt.dPol.firstRepayInd[j,i]=1
            sExt.dPol.defThreshold[j,i]=s.income.mBounds[1]

            sExt.polExtension.mListLength[j,i]=1
        else
            sExt.polExtension.mAPThreshold[i][1,j]=s.income.mBounds[2]
            sExt.polExtension.apPolicy[i][1,j]=lb
            #When the government sometimes defaults in income/debt state (i,j), further set:
            #3. the first index of the borrowing policy function which is relevant for the integration to 1
            #4. the threshold for default to uInverse(initial value of default - discounted continuation value of repayment under the optimal level of next period borrowing) minus
            #the level of consumption when m=0 under repayment at the optimal level of next period borrowing
            #5. the length of the list of thresholds for this state to 1
            sExt.dPol.firstRepayInd[j,i]=1
            sExt.dPol.defThreshold[j,i]=uInverse(m,sExt.VFD.vDInitGrid[i]-m.beta*sExt.VF.EVGrid[lb,i])-c(m,s,sExt,i,j,lb)

            sExt.polExtension.mListLength[j,i]=1
        end
    else
        #When lb<ub, the borrowing policy function is more complicated. In this case, we initialize three temporary variables to be used in a while loop. They are:
        #1. the length of the borrowing policy function should the while loop terminate
        #2. the index of the level of borrowing which at position tempLength in the borrowing policy function for which we want to find the threshold level of m below
        #which it is chosen
        #3. the index of the level of borrowing which is chosen above the threshold level of m in 2.
        #We begin by setting the length to one, the index of level of borrowing associated with that position to index chosen at the lowest possible value
        #of the m shock, and the alternative level of borrowing to 1 less than that index
        tempLength=1
        tempOldAP=lb
        tempNewAP=lb+1
        #Note that, in the range which is used for this function ub is the unrestricted maximizer when m=m_min and lb is the unrestricted maximizer when m=m_max
        #While we are still considering levels of borrowing which fall in the relevant range, continue the loop.
        while tempNewAP<=ub

            #Get the level of consumption at the current two levels of borrowing under consideration when m=0
            tempCM0Old=c(m,s,sExt,i,j,tempOldAP)
            tempCM0New=c(m,s,sExt,i,j,tempNewAP)

            #If borrowing less results in higher consumption, immediately "remove" the "old" value from the policy function and deincrement the length variable.
            #Otherwise, proceed to find the threshold.

            #A second condition here has been added due to an edge case that arose when finding the indifference levels
            #of beta for a consumer. In general, if the bond price function has been calculated correctly, asset choice
            #should be increasing in current assets (for a fixed value of income). This assumes, however,
            #arbitrary precision of all values involved, which we do not in general
            #use when programming. Roundoff errors resulted in cases where the second condition held as true in only
            #the first iteration after specifying a new value of beta for the government. Since, if we do not include
            #the second condition (which is rather strict; note the equality requirement for consumption), tempLength-1
            #can be zero, the update to tempOldAP will throw and error. To change this, uncomment the if statement just above the specification
            #of tempEVDiff, comment out the next two lines after the following comment starting with "Calculate",
            #and uncomment the specification of tempEVDiff which occurs in the following else block.

            #if tempCM0New>=tempCM0Old
            #Calculate the term (beta*(V_Old-V_New))
            tempEVDiff=m.beta*(sExt.VF.EVGrid[tempOldAP,i]-sExt.VF.EVGrid[tempNewAP,i])
            if (tempCM0New>tempCM0Old)||((tempCM0New==tempCM0Old)&&(tempEVDiff<0.0))
                #if tempLength==1
                #    println([j,i])
                #    println([tempCM0Old,tempCM0New])
                #    println(tempEVDiff)
                #end
                if tempLength!=1
                    tempOldAP=sExt.polExtension.apPolicy[i][tempLength-1,j]
                    tempLength+=-1
                else
                    sExt.dPol.attemptedSqRootNeg[j,i]=true
                    return sExt
                end


            else
                #Here we solve the equation 1/(C_New+m)-1/(C_Old+m)+beta*(V_Old-V_New)=0 for m. Some algebra yields:
                #(C_Old-C_New)/((C_New+m)*(C_Old+m))+beta*(V_Old-V_New)=0
                #(C_New+m)*(C_Old+m)+(C_Old-C_New)/(beta*(V_Old-V_New)=0
                #m^2+(CNew+C_Old)*m+C_New*C_Old+(C_Old-C_New)/(beta*(V_Old-V_New))=0
                #This is a parabola which opens up. We are interested in the rightmost root (the other one is, I believe, in general below whatever previous threshold was found)
                #We use the quadratic formula to do so below.

                #Calculate the term (beta*(V_Old-V_New))
                #tempEVDiff=m.beta*(sExt.VF.EVGrid[tempOldAP,i]-sExt.VF.EVGrid[tempNewAP,i])

                #tempEVDiff should, in general always be negative (if Z(y,a') is increasing in a'). If it is 0.0, set it to a negative value of very small magnitude.
                if tempEVDiff==0.0
                    tempEVDiff=-eps(F)
                end
                #if ((tempCM0New-tempCM0Old)/tempEVDiff-tempCM0New*tempCM0Old+0.25*(tempCM0New+tempCM0Old)^2)<0.0
                #    println([j,i])
                #    println([lb,ub])
                #    println([tempOldAP,tempNewAP])
                #    println([tempCM0Old,tempCM0New])
                #    println([sExt.VF.EVGrid[tempOldAP,i],sExt.VF.EVGrid[tempNewAP,i]])
                #    println([tempEVDiff])
                #    println([(tempCM0New-tempCM0Old)/tempEVDiff-tempCM0New*tempCM0Old+0.25*(tempCM0New+tempCM0Old)^2])
                #    println([(tempCM0New-tempCM0Old)/tempEVDiff,tempCM0New*tempCM0Old,0.25*(tempCM0New+tempCM0Old)^2])
                #end
                if ((tempCM0New-tempCM0Old)/tempEVDiff-tempCM0New*tempCM0Old+0.25*(tempCM0New+tempCM0Old)^2)<0.0
                    sExt.dPol.attemptedSqRootNeg[j,i]=true
                    return sExt
                end

                #Calculate the threshold value of m at which the government is indifferent between to the two levels of borrowing
                tempMThres=-0.5*(tempCM0New+tempCM0Old)+sqrt((tempCM0New-tempCM0Old)/tempEVDiff-tempCM0New*tempCM0Old+0.25*(tempCM0New+tempCM0Old)^2)

                #If tempLength is 1 (and therefore tempOldAP=lb; this will be the case even when we are not in the first iteration of the while loop):
                #1. Set the first entry of the threshold component of the borrowing policy function appropriately
                #2. Set the first entry of the borrowing policy function indices appropriately
                #3. Set the index for which we want to determine the upper bound of m values such that it is chosen to the index of the value for which we just found the threshold
                #4. Set the index of the borrowing value which should be used as the alternative in the indifference calculation for the next threshold to be one less than the index in 3.
                #5. Increment the length counter
                if tempLength==1
                    sExt.polExtension.mAPThreshold[i][tempLength,j]=tempMThres
                    sExt.polExtension.apPolicy[i][tempLength,j]=tempOldAP
                    tempOldAP=tempNewAP
                    tempNewAP=tempOldAP+1
                    tempLength+=1
                else
                    #If the borrowing policy function currently contains more entries than just what occurs at the lowest realization of m, we first check whether the threshold just
                    #calculated is less than the threshold above which tempOldAP SHOULD be chosen. If that is the case, then tempOldAP is NEVER chosen. We set tempOldAP to the value of
                    #the entry just above it in the borrowing policy function and deincrement the length counter (so that tempOldAP's presence in the borrowing policy function will
                    #eventually be overwritten and replaced with the proper value). Note that in this case tempNewAP is not changed.
                    if tempMThres<sExt.polExtension.mAPThreshold[i][tempLength-1,j]
                        tempOldAP=sExt.polExtension.apPolicy[i][tempLength-1,j]
                        tempLength+=-1
                    else
                        #If the new threshold is weakly greater than the old one:
                        #1. we add the upper bound of m values at which tempOldAP is chosen to the borrowing policy function
                        #2. record that it is chosen in the relevant interval
                        #3. make the new index for which an upper bound is to be determined the one which was used to find the upper bound in 1.
                        #4. set the index which we will attempt to used to find an upper bound for the m values at which the index in #3. is chosen to 1 less than the index in #3.
                        #5. increment the length counter.
                        sExt.polExtension.mAPThreshold[i][tempLength,j]=tempMThres
                        sExt.polExtension.apPolicy[i][tempLength,j]=tempOldAP
                        tempOldAP=tempNewAP
                        tempNewAP=tempOldAP+1
                        tempLength+=1
                    end
                end
            end
        end
        #When the above while loop exits, it will always be the case that tempOldAP=lb (since lb was the unrestricted maximizer at m=m_max) and tempNewAP
        #was strictly less than lb. We know the upper bound of values at which lb is chosen, and we know that it IS chosen. Record those facts, and record that
        #the first tempLength indices of the borrowing policy function are relevant for the current iteration's solution
        sExt.polExtension.mAPThreshold[i][tempLength,j]=s.income.mBounds[2]
        sExt.polExtension.apPolicy[i][tempLength,j]=ub
        sExt.polExtension.mListLength[j,i]=tempLength


        #We now move on determining default thresholds, if necessary
        if sExt.dPol.neverDefault[j,i]==true
            #If the government never defaults:
            #1. set the default threshold to the lower bound of the m shock
            #2. set the value of the first index of the borrowing policy function relevant for the integration steps to 1
            #3. set the probability of default to 0
            #4. calculate the probabilities associated with each choice of next period borrowing
            sExt.dPol.defThreshold[j,i]=s.income.mBounds[1]
            sExt.dPol.firstRepayInd[j,i]=1

        else
            #When the government does default for some values of m but not for others, determine which of the intervals specified by s.income.mBounds and s.pol.mAPThreshold it defaults in
            #Initialize this index at one
            defIdx=1
            #Continue until we reach the last relevant index
            #At each step, check whether the value at the upper bound of m values associated with that interval is weakly greater than the value of default. Once this is true, break the loop.
            #If it is false, increment the index of the interval in which we need to check whether default occurs
            #THESE STEPS ALL ASSUME THAT THIS SECTION IS NEVER CALLED WHEN EITHER alwaysDefault[j,i] OR neverDefault[j,i] ARE TRUE.
            while defIdx<=sExt.polExtension.mListLength[j,i]
                tempC=c(m,s,sExt,i,j,sExt.polExtension.apPolicy[i][defIdx,j],sExt.polExtension.mAPThreshold[i][defIdx,j])
                if u(m,tempC)+m.beta*sExt.VF.EVGrid[sExt.polExtension.apPolicy[i][defIdx,j],i]>=sExt.VFD.vDInitGrid[i]
                    break
                else
                    defIdx+=1
                end
            end

            #Set the index of the first element of the borrowing policy function relevant for the integration steps to the index of the interval in which default occurs
            sExt.dPol.firstRepayInd[j,i]=defIdx

            #Calculate the threshold at which default occurs within that interval
            sExt.dPol.defThreshold[j,i]=uInverse(m,sExt.VFD.vDInitGrid[i]-m.beta*sExt.VF.EVGrid[sExt.polExtension.apPolicy[i][defIdx,j],i])-c(m,s,sExt,i,j,sExt.polExtension.apPolicy[i][defIdx,j])

        end
    end

    return sExt
end


#solveRepayRow! solves the government's problem under repayment for every possible incoming level of borrowing
#given a specific value for the persistent component of income. It takes as its arguments:
#1. m: a model specification
#2. s: the collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
#4. i: the index of the value of the persistent component of income
#5. setExactP: a Boolean variable indicating whether or not to calculate exact transition probabilities

#Its output is just modified versions of s and g2

function solveRepayRow!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S},i::S) where{F<:Real,S<:Integer,T}

    tempAInd=1

    #initialize temporary variables for the upper and lower bounds of the grid search as well as the index of current asset holdings
    tempLB=1
    tempUB=sExt.pol.apLowM[sExt.aExtPoints+1,i]

    #Solve the government's problem at the lowest value of incoming borrowing under the lowest m shock
    sExt.polExtension.apLowM[1,i],tempVLow,g2.feasibleSolutionL[1,i]=solveRepayChoice(m,s,sExt,i,1,s.income.mBounds[1],tempLB,tempUB)

    #Mark that the problem has indeed been solved at the lowest value of incoming borrowing under the lowest m shock
    g2.solveMarkL[1,i]=true

    if g2.feasibleSolutionL[1,i]==true
        tempLB=sExt.polExtension.apLowM[1,i]
    end

    sExt.polExtension.apHighM[1,i],tempVHigh,g2.feasibleSolutionH[1,i]=solveRepayChoice(m,s,sExt,i,1,s.income.mBounds[2],tempLB,tempUB)

    #Mark that the problem has indeed been solved at the lowest value of incoming borrowing under the highest m shock
    g2.solveMarkH[1,i]=true

    #Set the alwaysDefault and neverDefault variables appropriately
    if tempVLow>=sExt.VFD.vDInitGrid[i]
        sExt.dPol.neverDefault[1,i]=true
        sExt.dPol.alwaysDefault[1,i]=false
    elseif tempVHigh<sExt.VFD.vDInitGrid[i]
        sExt.dPol.neverDefault[1,i]=false
        sExt.dPol.alwaysDefault[1,i]=true
    else
        sExt.dPol.neverDefault[1,i]=false
        sExt.dPol.alwaysDefault[1,i]=false
    end


    #reinitialize temporary variables for the upper and lower bounds of the grid search
    tempLB=1
    tempUB=sExt.pol.apHighM[sExt.aExtPoints+1,i]

    #Perform the exact same steps for the highest level of incoming borrowing
    sExt.polExtension.apLowM[sExt.aExtPoints,i],tempVLow,g2.feasibleSolutionL[sExt.aExtPoints,i]=solveRepayChoice(m,s,sExt,i,sExt.aExtPoints,s.income.mBounds[1],tempLB,tempUB)

    g2.solveMarkL[sExt.aExtPoints,i]=true


    if g2.feasibleSolutionL[sExt.aExtPoints,i]==true
        tempLB=sExt.polExtension.apLowM[sExt.aExtPoints,i]
    end

    sExt.polExtension.apHighM[sExt.aExtPoints,i],tempVHigh,g2.feasibleSolutionH[sExt.aExtPoints,i]=solveRepayChoice(m,s,sExt,i,sExt.aExtPoints,s.income.mBounds[2],tempLB,tempUB)

    g2.solveMarkH[sExt.aExtPoints,i]=true

    if tempVLow>=sExt.VFD.vDInitGrid[i]
        sExt.dPol.neverDefault[sExt.aExtPoints,i]=true
        sExt.dPol.alwaysDefault[sExt.aExtPoints,i]=false
    elseif tempVHigh<sExt.VFD.vDInitGrid[i]
        sExt.dPol.neverDefault[sExt.aExtPoints,i]=false
        sExt.dPol.alwaysDefault[sExt.aExtPoints,i]=true
    else
        sExt.dPol.neverDefault[sExt.aExtPoints,i]=false
        sExt.dPol.alwaysDefault[sExt.aExtPoints,i]=false
    end

    #Loop over the number of levels in the search parameters
    for levelInd in 1:(sExt.apSearchParams.numLevels)
        #Loop over the number of points in the current level
        for pointInd in 1:(sExt.apSearchParams.levelPoints[levelInd])
            #Set the index of the point currently under consideration
            tempAInd=sExt.apSearchParams.pointGrids[levelInd][pointInd]
            #If the problem has not been solved at this index for the highest realization of m, solve it. Otherwise (i.e. if the government was found to have defaulted at a
            #lower debt level under the highest realization of m), skip solving it, since the result does not matter.
            if g2.solveMarkH[tempAInd,i]==false

                #Get the indices associated with points for which the problem has already been solved whose solutions, if feasible, form upper and lower bounds
                #for the solution associated with the current index.
                tempUBInd=sExt.apSearchParams.pointUB[levelInd][pointInd]
                tempLBInd=sExt.apSearchParams.pointLB[levelInd][pointInd]

                #Check whether a feasible solution was found at those two points. If so, use the best bound available. Otherwise, use the alternative bounds already known.
                tempUB=ifelse(g2.feasibleSolutionH[tempUBInd,i]==true,sExt.polExtension.apHighM[tempUBInd,i],sExt.pol.apHighM[sExt.aExtPoints+1,i])
                tempLB=ifelse(g2.feasibleSolutionH[tempLBInd,i]==true,sExt.polExtension.apHighM[tempLBInd,i],1)

                #Solve the problem, check whether the solution resulted in strictly positive consumption, and mark that we have solved the problem at the current index
                sExt.polExtension.apHighM[tempAInd,i],tempVHigh,g2.feasibleSolutionH[tempAInd,i]=solveRepayChoice(m,s,sExt,i,tempAInd,s.income.mBounds[2],tempLB,tempUB)

                g2.solveMarkH[tempAInd,i]=true

                #If the value to the government under the highest realization of the m shock is strictly less than the initial value of default, then this will be true for:
                #1. all lower values of the m shock at this level of income and incoming asset.
                #2. all higher values of incoming asset at this level of income and the highest m shock.
                #3. all higher values of incoming asset at any level of the m shock (applying 1 to the states in 2 yields this result)
                #We use this information to mark that we know government default policy at such points (and because the government always defaults do not need
                #to solve the government's problem under repayment for any value of the m shock.
                #We also update the minAlwaysDefInd variable so that, should this condition be true for another value of debt, we need not update values which have already been set correctly.

                #If the value to the government under the highest realization of the m shock is weakly greater than the initial value of default, then we know that there are at least some values of
                #m at which the government does not default, so we set the always default variable appropriately.
                if tempVHigh<sExt.VFD.vDInitGrid[i]
                    sExt.dPol.neverDefault[g2.maxAlwaysDefInd[i]:tempAInd,i].=false
                    sExt.dPol.alwaysDefault[g2.maxAlwaysDefInd[i]:tempAInd,i].=true
                    g2.solveMarkH[g2.maxAlwaysDefInd[i]:tempAInd,i].=true
                    g2.solveMarkL[g2.maxAlwaysDefInd[i]:tempAInd,i].=true
                    g2.maxAlwaysDefInd[i]=max(g2.maxAlwaysDefInd[i],tempAInd)
                else
                    sExt.dPol.alwaysDefault[tempAInd,i]=false
                end
            end
            #If the problem has not been solved at this index for the lowest realization of m, solve it. Otherwise (i.e. if the government was found to have defaulted at a
            #lower debt level under the highest realization of m), skip solving it, since the result does not matter.
            if g2.solveMarkL[tempAInd,i]==false
                #Get the indices associated with points for which the problem has already been solved whose solutions, if feasible, form upper and lower bounds
                #for the solution associated with the current index.
                tempUBInd=sExt.apSearchParams.pointUB[levelInd][pointInd]
                tempLBInd=sExt.apSearchParams.pointLB[levelInd][pointInd]

                #Check whether a feasible solution was found at those two points. If so, use the best bound available. Otherwise, use the alternative bounds already known.
                tempUB=ifelse(g2.feasibleSolutionL[tempUBInd,i]==true,sExt.polExtension.apLowM[tempUBInd,i],sExt.pol.apLowM[sExt.aExtPoints+1,i])
                tempLB=ifelse(g2.feasibleSolutionL[tempLBInd,i]==true,sExt.polExtension.apLowM[tempLBInd,i],1)

                #Solve the problem, check whether the solution resulted in strictly positive consumption, and mark that we have solved the problem at the current index
                sExt.polExtension.apLowM[tempAInd,i],tempVLow,g2.feasibleSolutionL[tempAInd,i]=solveRepayChoice(m,s,sExt,i,tempAInd,s.income.mBounds[1],tempLB,tempUB)

                g2.solveMarkL[tempAInd,i]=true

                #Check whether the value to the government under the lowest possible value of the m shock is strictly greater than the initial value of default.
                #If this is the case, then the government never defaults in this state, and we can set the neverDefault variable to true. Otherwise, there are at least some values
                #of m at which the government defaults, so we set it to false.
                if tempVLow>=sExt.VFD.vDInitGrid[i]
                    sExt.dPol.neverDefault[tempAInd,i]=true
                else
                    sExt.dPol.neverDefault[tempAInd,i]=false
                end
            end
        end
    end

    #Loop over levels of debt.
    for j in 1:sExt.aExtPoints
        #Whenever the government does not always default (so we care about its full policy functions), construct the full policy function and probability objects associated with that state.
        if sExt.dPol.alwaysDefault[j,i]==false
            solveRepayInterior!(m,s,sExt,i,j)
        end
    end


    return s,g2

end


#integrateMVApprox performs the integration step in V(y,b)=E_m[W(y,m,b)] for a specific value of the persistent component of income and incoming asset, given government policies.
#This function performs the integration exactly as described in Chatterjee and Eyigungor (2012). Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. i: the index of the persistent component of income
#4. j: the index of the incoming level of borrowing
#It returns a single number which is the value of the integral.

#The structure and internal logic of this function are essentially identical to those of the one directly above it.
function integrateMVAlmostExact(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},i::S,j::S) where{F<:Real,S<:Integer,T}

    #The structure and internal logic of this function are essentially identical to those of the one directly above it.
    outValV=0.0
    if sExt.dPol.alwaysDefault[j,i]==true
        outValV+=sExt.VFD.vDInitGrid[i]
    elseif sExt.dPol.neverDefault[j,i]==true

        #If the government never defaults, then we begin the integration over the result under repayment at the first interval.
        #We set the location of the current upper bound in the m grid to the second point, the index of the relevant borrowing policy entry to the first,
        #and the most recent value of m to its lower bound
        mUBInd=2
        aPolicyInd=1
        lastMVal=s.income.mGrid[1]
        tempCBase=0.0
        tempCLow=0.0
        tempCHigh=0.0
        #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
        #it should never matter, and if it does, the process should in general not lead to convergence).
        while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(sExt.polExtension.mListLength[j,i]))
            #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
            #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
            #policy function
            if sExt.polExtension.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                tempCBase=c(m,s,sExt,i,j,sExt.polExtension.apPolicy[i][aPolicyInd,j])
                tempCLow=tempCBase+lastMVal
                tempCHigh=tempCBase+sExt.polExtension.mAPThreshold[i][aPolicyInd,j]
                outValV+=s.income.mProb[mUBInd-1]*((log(tempCLow)-log(tempCHigh))/s.income.mRes+(sExt.polExtension.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*m.beta*sExt.VF.EVGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i])
                lastMVal=sExt.polExtension.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                tempCBase=c(m,s,sExt,i,j,sExt.polExtension.apPolicy[i][aPolicyInd,j])
                tempCLow=tempCBase+lastMVal
                tempCHigh=tempCBase+s.income.mGrid[mUBInd]
                #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                #and then increment only the variable marking our position in the m grid
                outValV+=s.income.mProb[mUBInd-1]*((log(tempCLow)-log(tempCHigh))/s.income.mRes+(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*m.beta*sExt.VF.EVGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i])
                lastMVal=s.income.mGrid[mUBInd]

                mUBInd+=1
            end
        end

    else
        #If the government does default only sometimes, find the position in the m grid of the upper bound of the interval in which the threshold level of m for default falls
        #and set the last value of m observed to the value directly preceding that upper bound
        mUBInd=searchsortedfirst(s.income.mGrid,sExt.dPol.defThreshold[j,i])
        lastMVal=s.income.mGrid[mUBInd-1]
        #If there are any entire intervals in which default occurs, add their contribution to the output value
        if mUBInd>2
            for k in 3:mUBInd
                outValV+=s.income.mProb[k-2]*sExt.VFD.vDInitGrid[i]

            end
        end
        #Set the index of the first relevant entry of the borowing policy function
        aPolicyInd=sExt.dPol.firstRepayInd[j,i]

        #Add the contribution of default in the interval in which the threshold lies to the output value
        outValV+=s.income.mProb[mUBInd-1]*(sExt.dPol.defThreshold[j,i]-lastMVal)/s.income.mRes*sExt.VFD.vDInitGrid[i]

        #Update the last value of m observed to the threshold level of m at which default occurs
        lastMVal=sExt.dPol.defThreshold[j,i]
        tempCBase=0.0
        tempCLow=0.0
        tempCHigh=0.0
        #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
        #it should never matter, and if it does, the process should in general not lead to convergence).
        while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(sExt.polExtension.mListLength[j,i]))
            #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
            #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
            #policy function
            if sExt.polExtension.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                tempCBase=c(m,s,sExt,i,j,sExt.polExtension.apPolicy[i][aPolicyInd,j])
                tempCLow=tempCBase+lastMVal
                tempCHigh=tempCBase+sExt.polExtension.mAPThreshold[i][aPolicyInd,j]
                outValV+=s.income.mProb[mUBInd-1]*((log(tempCLow)-log(tempCHigh))/s.income.mRes+(sExt.polExtension.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*m.beta*sExt.VF.EVGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i])
                lastMVal=sExt.polExtension.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                tempCBase=c(m,s,sExt,i,j,sExt.polExtension.apPolicy[i][aPolicyInd,j])
                tempCLow=tempCBase+lastMVal
                tempCHigh=tempCBase+s.income.mGrid[mUBInd]
                #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                #and then increment only the variable marking our position in the m grid
                outValV+=s.income.mProb[mUBInd-1]*((log(tempCLow)-log(tempCHigh))/s.income.mRes+(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*m.beta*sExt.VF.EVGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i])
                lastMVal=s.income.mGrid[mUBInd]
                mUBInd+=1
            end
        end
    end

    return outValV
end


#updateV! performs the full integration across states to calculate V(y,b)=E_m[W(y,m,b)]. Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
#4. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of vNew
function updateV!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S})  where{F<:Real,S<:Integer,T}

    Threads.@threads for i in 1:(m.yParams.yPoints-1)
    #for i in 1:(m.yParams.yPoints-1)
        for j in 1:sExt.aExtPoints
            g2.VFExtension.vGrid[j,i]=integrateMVAlmostExact(m,s,sExt,i,j)
        end
    end
    return g2
end


#updateEV! performs the second integration step in updating the main value function. Specifically, it calculates Z(y,b')=E_y'[V(y',b')].
#Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
function updateEV!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S})  where{F<:Real,S<:Integer,T}
    mul!(g2.VFExtension.EVGrid,g2.VFExtension.vGrid,s.income.yTMat)

    return g2
end


#integrateMQSumApprox performs the first step of the integration in the functional equation which q must satisfy for a single value of the persistent component of income
#and incoming borrowing level. It returns the term E_m[(1-d(y,m,b))*(lambda+(1-lambda)(coup+q(y,m,b'(y,m,b))))]. This implementation uses exactly the approximation
#described in Chatterjee and Eyigungor (2012). Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. i: the index of the persistent component of income
#4. j: the index of the incoming level of borrowing
#It returns a single number which is the value of the integral.

#The structure and internal logic of this function is essentially identical to that of integrateMVApprox
function integrateMQSumApprox(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},i::S,j::S) where{F<:Real,S<:Integer,T}
    #Initialize the output value
    outValQ=0.0
    #tempPSum=0.0
    #If the government ever repays in the given state, sum the appropriate values
    if sExt.dPol.neverDefault[j,i]==true
        mUBInd=2
        aPolicyInd=1
        lastMVal=s.income.mGrid[1]
        while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(sExt.polExtension.mListLength[j,i]))
            if (sExt.polExtension.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                outValQ+=s.income.mProb[mUBInd-1]*(sExt.polExtension.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(m.lambda+(1.0-m.lambda)*(m.coup+sExt.qGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i]))
                #tempPSum+=s.income.mProb[mUBInd-1]*(sExt.polExtension.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                lastMVal=sExt.polExtension.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                outValQ+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(m.lambda+(1.0-m.lambda)*(m.coup+sExt.qGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i]))
                #tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                lastMVal=s.income.mGrid[mUBInd]
                mUBInd+=1
            end
        end
        #outValQ*=tempPSum^(-1)
    elseif sExt.dPol.alwaysDefault[j,i]==false

        aPolicyInd=sExt.dPol.firstRepayInd[j,i]
        lastMVal=sExt.dPol.defThreshold[j,i]

        mUBInd=searchsortedfirst(s.income.mGrid,lastMVal)

        #if mUBInd>2
        #    for k in 3:mUBInd
        #        tempPSum+=s.income.mProb[k-2]
        #    end
        #end
        #tempPSum+=s.income.mProb[mUBInd-1]*(sExt.dPol.defThreshold[j,i]-s.income.mGrid[mUBInd-1])/s.income.mRes

        while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(sExt.polExtension.mListLength[j,i]))
            if (sExt.polExtension.mAPThreshold[i][aPolicyInd,j])<(s.income.mGrid[mUBInd])
                outValQ+=s.income.mProb[mUBInd-1]*(sExt.polExtension.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*(m.lambda+(1.0-m.lambda)*(m.coup+sExt.qGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i]))
                #tempPSum+=s.income.mProb[mUBInd-1]*(sExt.polExtension.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes
                lastMVal=sExt.polExtension.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                outValQ+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*(m.lambda+(1.0-m.lambda)*(m.coup+sExt.qGrid[sExt.polExtension.apPolicy[i][aPolicyInd,j],i]))
                #tempPSum+=s.income.mProb[mUBInd-1]*(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes
                lastMVal=s.income.mGrid[mUBInd]
                mUBInd+=1
            end
        end
        #outValQ*=tempPSum^(-1)
    end

    return outValQ

end



#updateQ! performs the full integration across states to calculate q(y,b'). Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
#4. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of qNew
function updateQ!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S}) where{F<:Real,S<:Integer,T}

    #Perform the first step of the integration
    Threads.@threads for i in 1:(m.yParams.yPoints)
    #for i in 1:(m.yParams.yPoints)
        for j in 1:(sExt.aExtPoints)
            g2.qSumExtension[j,i]=integrateMQSumApprox(m,s,sExt,i,j)/m.R
        end
    end


    #Perform the second step of the integration, and store the result in qNew before returning it
    mul!(g2.qGridExtension,g2.qSumExtension,s.income.yTMat)

    #copyto!(g2.qGrid,(g2.qSum*s.income.yTMat))
    return g2
end



#vfiGOneStep! is the main function of the for solving method for the model in Chatterjee and Eyigungor (2012)
#1. m: a model specification
#2. s: a collection of objects used to solve the model which will be modified
#3. tol: the tolerance for convergence for both the value functions and the bond price function
#4. maxIter: the maximum number of iterations
#5. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of s
function vfiGOneStep!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},tol::F,maxIter::S,mixFacQ::F,mixFacV::F) where{F<:Real,S<:Integer,T}




    #Initialize new copies of the value functions and the bond price function
    g2=makeLTBExtensionUpdate(m,s,sExt)


    #Initialize  a counter for the number of iterations and measures of sup norm distance between various objects
    iCount=1
    vFDist=tol+1.0
    EVDist=tol+1.0
    qDist=tol+1.0

    vDist=tol+1.0
    maxDist=tol+1.0

    somethingIsNonMonotone=false

    #Iterate until the maximum number of iterations has been reached or the sup norm distance between successive iterations drops below the tolerance for convergence
    while (iCount<=maxIter)&(maxDist>=tol)

        #Reset all the objects used to speed up the grid search step
        g2.solveMarkL.=false
        g2.solveMarkH.=false
        g2.feasibleSolutionL.=false
        g2.feasibleSolutionH.=false
        g2.maxAlwaysDefInd.=1
        g2.maxAlwaysDefInd[m.yParams.yPoints]=sExt.aExtPoints


        #Iterate over states of the persistent component of income
        Threads.@threads for i in 1:(m.yParams.yPoints-1)
        #for i in 1:(m.yParams.yPoints-1)
            #solve the problem at the current state of the persistent component of income
            solveRepayRow!(m,s,sExt,g2,i)
        end

        somethingIsNonMonotone=maximum(sExt.dPol.attemptedSqRootNeg)

        if somethingIsNonMonotone==true
            break
        end

        #Perform the integration steps to obtain new versions of V, Z, and q
        updateV!(m,s,sExt,g2)
        updateEV!(m,s,sExt,g2)
        updateQ!(m,s,sExt,g2)


        #Calculate the measures of distance between successive iterations for V, Z, and q
        vFDist=maximum(abs.(g2.VFExtension.vGrid.-sExt.VFExtension.vGrid))
        EVDist=maximum(abs.(g2.VFExtension.EVGrid.-sExt.VFExtension.EVGrid))
        qDist=maximum(abs.(g2.qGridExtension.-sExt.qGrid[1:sExt.aExtPoints,:]))

        #Update the guesses for V, Z, and q
        sExt.VFExtension.vGrid.=mixFacV*sExt.VFExtension.vGrid.+(1.0-mixFacV)*g2.VFExtension.vGrid
        sExt.VFExtension.EVGrid.=mixFacV*sExt.VFExtension.EVGrid.+(1.0-mixFacV)*g2.VFExtension.EVGrid

        sExt.VF.vGrid[1:sExt.aExtPoints,:].=sExt.VFExtension.vGrid
        sExt.VF.EVGrid[1:sExt.aExtPoints,:].=sExt.VFExtension.EVGrid

        sExt.qGrid[1:sExt.aExtPoints,:].=mixFacQ*sExt.qGrid[1:sExt.aExtPoints,:].+(1.0-mixFacQ)*g2.qGridExtension

        copyto!(sExt.qSumExtension,g2.qSumExtension)

        #Calculate the maximum of the various distance measures for the value functions
        vDist=max(vFDist,EVDist)

        #Calculate the maximum of all the various distance measures
        maxDist=max(vDist,qDist)

        #At each iteration which is a multiple of 5% of the maximum number of iterations, print the current iteration number and variables tracking the convergence criteria
        if mod(iCount,100)==0
            # println([iCount,qDist,EVDist,vFDist])
            print(".")
        end



        #Increment the iteration counter
        iCount+=1


    end



    #Print final set of convergence statistics
    #println("\n", [iCount-1,qDist,EVDist,vFDist])

    #Update the marker tracking convergence of the value functions
    if vDist>tol
        sExt.noConvergeMark[1]=true
    else
        sExt.noConvergeMark[1]=false
    end

    for i in 1:(m.yParams.yPoints)
        for j in 1:(m.aPoints+sExt.aExtPoints)

            view(sExt.consM0Grid,:,j,i).=sExt.netRevM0A0Grid[j,i].-view(sExt.qGrid,:,i).*view(sExt.aGridIncr,:,j)

        end
    end
    sExt.pol.apLowM[1:sExt.aExtPoints,:].=sExt.polExtension.apLowM
    sExt.pol.apHighM[1:sExt.aExtPoints,:].=sExt.polExtension.apHighM
    sExt.pol.mListLength[1:sExt.aExtPoints,:].=sExt.polExtension.mListLength

    for i in 1:(m.yParams.yPoints)
        for j in 1:(sExt.aExtPoints)
            for jj in 1:(sExt.pol.mListLength[j,i])
                sExt.pol.mAPThreshold[i][jj,j]=sExt.polExtension.mAPThreshold[i][jj,j]
                sExt.pol.apPolicy[i][jj,j]=sExt.polExtension.apPolicy[i][jj,j]
            end
        end
    end
    #Return the modified collection of objects used to solve the model
    return sExt,g2,somethingIsNonMonotone,maxDist
end


#vfiGOneStep! is the main function of the for solving method for the model in Chatterjee and Eyigungor (2012)
#1. m: a model specification
#2. s: a collection of objects used to solve the model which will be modified
#3. tol: the tolerance for convergence for both the value functions and the bond price function
#4. maxIter: the maximum number of iterations
#5. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of s
function vfiGOneStepRM!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},tol::F,maxIter::S,mixFacQ::F,mixFacV::F) where{F<:Real,S<:Integer,T}




    #Initialize new copies of the value functions and the bond price function
    g2=makeLTBExtensionUpdate(m,s,sExt)


    #Initialize  a counter for the number of iterations and measures of sup norm distance between various objects
    iCount=1
    vFDist=tol+1.0
    EVDist=tol+1.0
    qDist=tol+1.0

    vDist=tol+1.0
    maxDist=tol+1.0

    somethingIsNonMonotone=false

    randMix=0.5
    mixFacQR=0.5
    mixFacVR=0.5

    #Iterate until the maximum number of iterations has been reached or the sup norm distance between successive iterations drops below the tolerance for convergence
    while (iCount<=maxIter)&(maxDist>=tol)

        #Reset all the objects used to speed up the grid search step
        g2.solveMarkL.=false
        g2.solveMarkH.=false
        g2.feasibleSolutionL.=false
        g2.feasibleSolutionH.=false
        g2.maxAlwaysDefInd.=1
        g2.maxAlwaysDefInd[m.yParams.yPoints]=sExt.aExtPoints


        #Iterate over states of the persistent component of income
        Threads.@threads for i in 1:(m.yParams.yPoints-1)
        #for i in 1:(m.yParams.yPoints-1)
            #solve the problem at the current state of the persistent component of income
            solveRepayRow!(m,s,sExt,g2,i)
        end

        somethingIsNonMonotone=maximum(sExt.dPol.attemptedSqRootNeg)

        if somethingIsNonMonotone==true
            break
        end

        #Perform the integration steps to obtain new versions of V, Z, and q
        updateV!(m,s,sExt,g2)
        updateEV!(m,s,sExt,g2)
        updateQ!(m,s,sExt,g2)

        randMix=rand()

        mixFacQR=randMix*mixFacQ+(1.0-randMix)*1.0
        mixFacVR=randMix*mixFacV+(1.0-randMix)*1.0

        #Calculate the measures of distance between successive iterations for V, Z, and q
        vFDist=maximum(abs.(g2.VFExtension.vGrid.-sExt.VFExtension.vGrid))
        EVDist=maximum(abs.(g2.VFExtension.EVGrid.-sExt.VFExtension.EVGrid))
        qDist=maximum(abs.(g2.qGridExtension.-sExt.qGrid[1:sExt.aExtPoints,:]))

        #Update the guesses for V, Z, and q
        sExt.VFExtension.vGrid.=mixFacVR*sExt.VFExtension.vGrid.+(1.0-mixFacVR)*g2.VFExtension.vGrid
        sExt.VFExtension.EVGrid.=mixFacVR*sExt.VFExtension.EVGrid.+(1.0-mixFacVR)*g2.VFExtension.EVGrid

        sExt.VF.vGrid[1:sExt.aExtPoints,:].=sExt.VFExtension.vGrid
        sExt.VF.EVGrid[1:sExt.aExtPoints,:].=sExt.VFExtension.EVGrid

        sExt.qGrid[1:sExt.aExtPoints,:].=mixFacQR*sExt.qGrid[1:sExt.aExtPoints,:].+(1.0-mixFacQR)*g2.qGridExtension

        copyto!(sExt.qSumExtension,g2.qSumExtension)

        #Calculate the maximum of the various distance measures for the value functions
        vDist=max(vFDist,EVDist)

        #Calculate the maximum of all the various distance measures
        maxDist=max(vDist,qDist)

        #At each iteration which is a multiple of 5% of the maximum number of iterations, print the current iteration number and variables tracking the convergence criteria
        if mod(iCount,100)==0
            # println([iCount,qDist,EVDist,vFDist])
            print(".")
        end



        #Increment the iteration counter
        iCount+=1


    end



    #Print final set of convergence statistics
    # println("\n", [iCount-1,qDist,EVDist,vFDist])

    #Update the marker tracking convergence of the value functions
    if vDist>tol
        sExt.noConvergeMark[1]=true
    else
        sExt.noConvergeMark[1]=false
    end

    for i in 1:(m.yParams.yPoints)
        for j in 1:(m.aPoints+sExt.aExtPoints)

            view(sExt.consM0Grid,:,j,i).=sExt.netRevM0A0Grid[j,i].-view(sExt.qGrid,:,i).*view(sExt.aGridIncr,:,j)

        end
    end
    sExt.pol.apLowM[1:sExt.aExtPoints,:].=sExt.polExtension.apLowM
    sExt.pol.apHighM[1:sExt.aExtPoints,:].=sExt.polExtension.apHighM
    sExt.pol.mListLength[1:sExt.aExtPoints,:].=sExt.polExtension.mListLength

    for i in 1:(m.yParams.yPoints)
        for j in 1:(sExt.aExtPoints)
            for jj in 1:(sExt.pol.mListLength[j,i])
                sExt.pol.mAPThreshold[i][jj,j]=sExt.polExtension.mAPThreshold[i][jj,j]
                sExt.pol.apPolicy[i][jj,j]=sExt.polExtension.apPolicy[i][jj,j]
            end
        end
    end
    #Return the modified collection of objects used to solve the model
    return sExt,g2,somethingIsNonMonotone,maxDist
end


#writeEssentials is a Quick function for saving the pieces of a model necessary for producing figures as well as the full set of objects
#which can be used to reconstruct the computed model given the specification only.
function writeEssentials(s::LongTermBondRDSaveEval{F,S},sExt::LongTermBondRDSaveExtension{F,S},filePrefix::String,fileDir::String) where{F<:Real,S<:Integer}
    writecsv(
        joinpath(fileDir, filePrefix*"_yGrid.csv"),
        s.income.yGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_yTMat.csv"),
        s.income.yTMat
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_yDefGrid.csv"),
        sExt.yDefGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_aGrid.csv"),
        sExt.aGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_qGrid.csv"),
        sExt.qGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_vGrid.csv"),
        sExt.VF.vGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_EVGrid.csv"),
        sExt.VF.EVGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_vDFutGrid.csv"),
        sExt.VFD.vDFutGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_vDInitGrid.csv"),
        sExt.VFD.vDInitGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_EVDGrid.csv"),
        sExt.VFD.EVDGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_vDFutFlow.csv"),
        sExt.VFD.vDFutFlow
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_vDInitFlow.csv"),
        sExt.VFD.vDInitFlow
    )
end


#writeEssentials is a Quick function for saving the pieces of a model necessary for producing figures as well as the full set of objects
#which can be used to reconstruct the computed model given the specification only.
function readModel!(m::LongTermBondRDSaveSpec{F,S},sExt::LongTermBondRDSaveExtension{F,S},filePrefix::String,fileDir::String) where{F<:Real,S<:Integer}
    sExt.qGrid.=readcsv(
        joinpath(fileDir, filePrefix*"_qGrid.csv")
    )
    sExt.VF.vGrid.=readcsv(
        joinpath(fileDir, filePrefix*"_vGrid.csv")
    )
    sExt.VF.EVGrid.=readcsv(
        joinpath(fileDir, filePrefix*"_EVGrid.csv")
    )

    sExt.VFExtension.vGrid.=sExt.VF.vGrid[1:sExt.aExtPoints,:]
    sExt.VFExtension.EVGrid.=sExt.VF.EVGrid[1:sExt.aExtPoints,:]
    for i in 1:(m.yParams.yPoints)
        for j in 1:(m.aPoints+sExt.aExtPoints)
            view(sExt.consM0Grid,:,j,i).=sExt.netRevM0A0Grid[j,i].-view(sExt.qGrid,:,i).*view(sExt.aGridIncr,:,j)
        end
    end
    return sExt
end


#updateQ! performs the full integration across states to calculate q(y,b'). Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
#4. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of qNew
function updateQSumRow!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S},i::S) where{F<:Real,S<:Integer,T}

    #Perform the first step of the integration

    for j in 1:(sExt.aExtPoints)
        g2.qSumExtension[j,i]=integrateMQSumApprox(m,s,sExt,i,j)/m.R
    end



    return g2
end


#Its output is simply a modified version of qNew
function updateQRow!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S},i::S) where{F<:Real,S<:Integer,T}

    #Perform the first step of the integration

    view(g2.qGridExtension,:,i).=sExt.qSumExtension*view(s.income.yTMat,:,i)



    return g2
end



#Its output is simply a modified version of vNew
function updateVRow!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S},i::S)  where{F<:Real,S<:Integer,T}


    for j in 1:sExt.aExtPoints

        g2.VFExtension.vGrid[j,i]=integrateMVAlmostExact(m,s,sExt,i,j)

    end



    return g2
end

#updateEV! performs the second integration step in updating the main value function. Specifically, it calculates Z(y,b')=E_y'[V(y',b')].
#Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
function updateEVRow!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},g2::LongTermBondRDExtensionUpdate{F,S},i::S)  where{F<:Real,S<:Integer,T}
    view(g2.VFExtension.EVGrid,:,i).=sExt.VFExtension.vGrid*view(s.income.yTMat,:,i)

    return g2
end



#Its output is simply a modified version of s
function vfiGOneStepRBR!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},sExt::LongTermBondRDSaveExtension{F,S},tol::F,maxIter::S,mixFacQ::F,mixFacQSum::F,mixFacEV::F,mixFacV::F) where{F<:Real,S<:Integer,T}




    #Initialize new copies of the value functions and the bond price function
    g2=makeLTBExtensionUpdate(m,s,sExt)


    #Initialize  a counter for the number of iterations and measures of sup norm distance between various objects
    iCount=1
    vFDist=tol+1.0
    EVDist=tol+1.0
    qDist=tol+1.0
    qSumDist=tol+1.0

    vDist=tol+1.0
    maxDist=tol+1.0

    #Iterate until the maximum number of iterations has been reached or the sup norm distance between successive iterations drops below the tolerance for convergence
    while (iCount<=maxIter)&(maxDist>=tol)

        #Reset all the objects used to speed up the grid search step
        g2.solveMarkL.=false
        g2.solveMarkH.=false
        g2.feasibleSolutionL.=false
        g2.feasibleSolutionH.=false
        g2.maxAlwaysDefInd.=1
        g2.maxAlwaysDefInd[m.yParams.yPoints]=sExt.aExtPoints

        qDist=0.0
        qSumDist=0.0
        vFDist=0.0
        EVDist=0.0

        #iRand=rand()
        #Iterate over states of the persistent component of income
        for iName in 1:(m.yParams.yPoints-1)
            i=ifelse(isodd(iCount),iName,m.yParams.yPoints-iName)
            #i=ifelse(iRand>0.5,iName,m.yParams.yPoints-iName)

            updateEVRow!(m,s,sExt,g2,i)

            EVDist=max(EVDist,maximum(abs.(view(g2.VFExtension.EVGrid,:,i).-view(sExt.VFExtension.EVGrid,:,i))))
            view(sExt.VFExtension.EVGrid,:,i).=mixFacEV.*view(sExt.VFExtension.EVGrid,:,i).+(1.0-mixFacEV).*view(g2.VFExtension.EVGrid,:,i)
            view(sExt.VF.EVGrid,1:sExt.aExtPoints,i).=view(sExt.VFExtension.EVGrid,:,i)

            updateQRow!(m,s,sExt,g2,i)

            qDist=max(qDist,maximum(abs.(view(g2.qGridExtension,:,i).-view(sExt.qGrid,1:sExt.aExtPoints,i))))
            view(sExt.qGrid,1:sExt.aExtPoints,i).=mixFacQ.*view(sExt.qGrid,1:sExt.aExtPoints,i).+(1.0-mixFacQ).*view(g2.qGridExtension,:,i)
            #view(sExt.qGrid,1:sExt.aExtPoints,i).=view(g2.qGridExtension,:,i)
            #solve the problem at the current state of the persistent component of income
            solveRepayRow!(m,s,sExt,g2,i)

            #Perform the integration steps to obtain new versions of V, Z, and q
            updateVRow!(m,s,sExt,g2,i)

            vFDist=max(vFDist,maximum(abs.(view(g2.VFExtension.vGrid,:,i).-view(sExt.VFExtension.vGrid,:,i))))
            view(sExt.VFExtension.vGrid,:,i).=mixFacV.*view(sExt.VFExtension.vGrid,:,i).+(1.0-mixFacV).*view(g2.VFExtension.vGrid,:,i)
            view(sExt.VF.vGrid,1:sExt.aExtPoints,i).=view(sExt.VFExtension.vGrid,:,i)

            updateQSumRow!(m,s,sExt,g2,i)

            qSumDist=max(qSumDist,maximum(abs.(view(g2.qSumExtension,:,i).-view(sExt.qSumExtension,:,i))))
            view(sExt.qSumExtension,:,i).=mixFacQSum.*view(sExt.qSumExtension,:,i).+(1.0-mixFacQSum).*view(g2.qSumExtension,:,i)


        end


        #Calculate the maximum of the various distance measures for the value functions
        vDist=max(vFDist,EVDist)

        #Calculate the maximum of all the various distance measures
        maxDist=max(vDist,qDist,qSumDist)

        #At each iteration which is a multiple of 5% of the maximum number of iterations, print the current iteration number and variables tracking the convergence criteria
        if mod(iCount,10)==0
            # println([iCount,qDist,qSumDist,EVDist,vFDist])
            print(".")
        end



        #Increment the iteration counter
        iCount+=1


    end



    #Print final set of convergence statistics
    # println("\n", [iCount-1,qDist,EVDist,vFDist])

    #Update the marker tracking convergence of the value functions
    if vDist>tol
        sExt.noConvergeMark[1]=true
    else
        sExt.noConvergeMark[1]=false
    end

    for i in 1:(m.yParams.yPoints)
        for j in 1:(m.aPoints+sExt.aExtPoints)

            view(sExt.consM0Grid,:,j,i).=sExt.netRevM0A0Grid[j,i].-view(sExt.qGrid,:,i).*view(sExt.aGridIncr,:,j)

        end
    end
    sExt.pol.apLowM[1:sExt.aExtPoints,:].=sExt.polExtension.apLowM
    sExt.pol.apHighM[1:sExt.aExtPoints,:].=sExt.polExtension.apHighM
    sExt.pol.mListLength[1:sExt.aExtPoints,:].=sExt.polExtension.mListLength

    for i in 1:(m.yParams.yPoints)
        for j in 1:(sExt.aExtPoints)
            for jj in 1:(sExt.pol.mListLength[j,i])
                sExt.pol.mAPThreshold[i][jj,j]=sExt.polExtension.mAPThreshold[i][jj,j]
                sExt.pol.apPolicy[i][jj,j]=sExt.polExtension.apPolicy[i][jj,j]
            end
        end
    end
    #Return the modified collection of objects used to solve the model
    return sExt,g2
end
