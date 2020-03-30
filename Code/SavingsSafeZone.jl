#The function makeIncomeProcess constructs a income process object based on the contents of a model specification.
#Its only argument is:
#1. m: a model specification

#Its output is simply a combinedIncomeProcessSave object

function makeIncomeProcess(m::LongTermBondRDSaveSpec{F,S}) where{F<:Real,S<:Integer}
    #Generate the grid of income values, the transition matrix, and the long run mean of the process
    yGrid,yTMat,yMean=tauchenRev(m.yParams)
    ydim=m.yParams.yPoints



    #Generate the distribution of m shock
    mDist=Distributions.Normal(m.mParams.mu,sqrt(m.mParams.epsilon2))

    #Set the bounds for the m shock
    mBounds=m.mParams.stdSpan*sqrt(m.mParams.epsilon2)*[-1.0,1.0].+m.mParams.mu

    #Create the grid of points for the m shock
    mGrid=collect(LinRange(mBounds[1],mBounds[2],m.mParams.mPoints))

    #Create the midPoints of the grids for the m shock
    mMidPoints=0.5*(mGrid[2:(m.mParams.mPoints)]+mGrid[1:(m.mParams.mPoints-1)])

    mRes=(mBounds[2]-mBounds[1])/(m.mParams.mPoints-1)

    #Create the probabilities for the intervals of the m shock
    mProb=zeros(m.mParams.mPoints-1)
    for i in 2:(m.mParams.mPoints)
        mProb[i-1]=cdf(mDist,mGrid[i])-cdf(mDist,mGrid[i-1])
    end
    tempMPSumInv=sum(mProb)^(-1)
    for i in 2:(m.mParams.mPoints)
        mProb[i-1]*=tempMPSumInv
    end


    #Calculate the inflation factor for the truncated distribution of m
    mInfl=(cdf(mDist,mBounds[2])-cdf(mDist,mBounds[1]))^(-1)


    return CombinedIncomeProcessSave(yGrid,yTMat,yMean,mBounds,mGrid,mMidPoints,mProb,mRes,mInfl,mDist)
end


#The function makeValueFunction constructs a value function object based on the contents of a model specification. Its only argument is:
#1. m: a model specification

#Its output is simply a vFuncSave object

function makeValueFunction(m::LongTermBondRDSaveSpec{F,S}) where{F<:Real,S<:Integer}
    #Set some aliases
    ydim=m.yParams.yPoints
    adim=m.aPoints
    #Initialize the value function arrays
    vGrid=zeros(adim,ydim)

    EVGrid=zeros(adim,ydim)





    return VFuncSave(vGrid,EVGrid)
end

#The function makePolicyGrids constructs a policy function object based on the contents of a model specification.
#Its only argument is:
#1. m: a model specification

#Its output is simply a policiesSave object

function makePolicyGrids(m::LongTermBondRDSaveSpec{F,S}) where{F<:Real,S<:Integer}
    #Set some aliases
    ydim=m.yParams.yPoints
    adim=m.aPoints
    #Initialize the collection of objects which define the government's policy functions
    apLowM=ones(S,adim,ydim)
    apHighM=ones(S,adim,ydim)
    mListLength=zeros(S,adim,ydim)
    mAPThreshold=Array{SparseMatrixCSC{F,S},1}(undef,ydim)
    apPolicy=Array{SparseMatrixCSC{S,S},1}(undef,ydim)


    sDummyF=sparse(zeros(F,adim,adim))
    sDummyS=sparse(zeros(S,adim,adim))
    for i in 1:ydim
        mAPThreshold[i]=deepcopy(sDummyF)
        apPolicy[i]=deepcopy(sDummyS)
    end

    return PoliciesSave(apLowM,apHighM,mListLength,mAPThreshold,apPolicy)
end

#The function makeLTBUpdate constructs a combinedIncomeProcessSave object based on the contents of a model specification and eval object.
#Its arguments are:
#1. m: a model specification
#2. s: the collection of objects used to solve the model

#Its output is simply a longTermBondRDSaveUpdate object

function makeLTBUpdate(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T}) where{F<:Real,S<:Integer,T}
    #Set some aliases
    ydim=m.yParams.yPoints
    adim=m.aPoints

    VF=deepcopy(s.VF)




    solveMarkH=zeros(Bool,adim,ydim)
    solveMarkL=zeros(Bool,adim,ydim)
    feasibleSolutionH=zeros(Bool,adim,ydim)
    feasibleSolutionL=zeros(Bool,adim,ydim)


    return LongTermBondRDSaveUpdate(VF,solveMarkH,solveMarkL,feasibleSolutionH,feasibleSolutionL)
end



#The function makeEval constructs a longTermBondRDSaveEval object based on the contents of a model specification. Its only argument is:
#1. m: a model specification

#Its output is simply:
#1. outputEval: a collection of objects to be used in solving the model

function makeEval(m::LongTermBondRDSaveSpec{F,S}) where{F<:Real,S<:Integer}

    #generate the income process, value function, and policy function objects
    income=makeIncomeProcess(m)
    VF=makeValueFunction(m)
    pol=makePolicyGrids(m)

    #Generate the grid of borrowing levels
    aGrid=collect(LinRange(m.aBounds[1],m.aBounds[2],m.aPoints))
    aGridIncr=aGrid.-(1.0-m.lambda).*aGrid'

    #Set some aliases
    ydim=m.yParams.yPoints
    adim=m.aPoints

    #Initialize the bond price function
    qBase=(m.lambda+(1-m.lambda)*m.coup)/(m.R-(1.0-m.lambda))



    #Find the location of 0.0 in the borrowing grid
    a0Ind=searchsortedfirst(aGrid,zero(F))

    #Check that the result is actually 0.0, i.e. the grid has not been misspecified as an evenly spaced grid which skips 0.0
    @assert aGrid[a0Ind]==zero(F)


    #q=1/R*(lambda+(1-lambda)*(coup+q)
    #(R-(1-lambda))*q=lambda+(1-lambda)*coup
    #q=(lambda+(1-lambda)*coup)/(R-(1-lambda))
    netRevM0A0Grid=zeros(adim,ydim)

    #Fill the grid of consumption values when m=0 and b'=(1-lambda)*b
    for i in 1:ydim
        for j in 1:adim
            netRevM0A0Grid[j,i]=income.yGrid[i]+(m.lambda+(1.0-m.lambda)*m.coup)*aGrid[j]
        end
    end

    #Initialize the grids for consumption when m=0 and b'=(1-lambda)*b and for when m=0 and b' takes a specified value
    consM0Grid=zeros(adim,adim,ydim)

    #Fill the grid of consumption values when m=0 and b'=(1-lambda)*b
    for i in 1:(m.yParams.yPoints)
        for j in 1:(m.aPoints)
            view(consM0Grid,:,j,i).=netRevM0A0Grid[j,i].-qBase.*view(aGridIncr,:,j)
        end
    end


    #Initialize the array of maximum feasible choices of savings
    maxAPInd=ones(S,adim,ydim)

    for i in 1:ydim
        for j in 1:adim
            apInd=m.aPoints
            maxAPInd[j,i]=adim
            while apInd>1
                if consM0Grid[apInd,j,i]+income.mBounds[2]>0.0
                    maxAPInd[j,i]=apInd
                    break
                else
                    apInd+=-1
                end
            end
        end
    end





    #Generate the set of obects used to speed the grid search for the highest and lowest values of b' chosen at every (y,b) state
    apSearchParams=genBisectSearchGrids(adim)

    #construct and return the longTermBondRDSaveEval object whose components are constructed above
    outputEval=LongTermBondRDSaveEval(income,VF,pol,aGrid,aGridIncr,qBase,netRevM0A0Grid,consM0Grid,maxAPInd,a0Ind,[true],apSearchParams)


    return outputEval

end




function makeEval(m::LongTermBondRDSaveSpec{F,S},aFineUB::F,aFinePoints::S) where{F<:Real,S<:Integer}

    @assert aFineUB<m.aBounds[2]
    @assert aFineUB>m.aBounds[1]
    @assert aFinePoints>0
    @assert aFinePoints<(m.aPoints-1)

    #generate the income process, value function, and policy function objects
    income=makeIncomeProcess(m)
    VF=makeValueFunction(m)
    pol=makePolicyGrids(m)

    aCoarsePoints=m.aPoints-aFinePoints
    aFineRes=(aFineUB-m.aBounds[1])/aFinePoints

    #Generate the grid of borrowing levels
    aGrid=vcat(collect(LinRange(m.aBounds[1],aFineUB-aFineRes,aFinePoints)),collect(LinRange(aFineUB,m.aBounds[2],aCoarsePoints)))
    aGridIncr=aGrid.-(1.0-m.lambda).*aGrid'

    #Set some aliases
    ydim=m.yParams.yPoints
    adim=m.aPoints

    #Initialize the bond price function
    qBase=(m.lambda+(1-m.lambda)*m.coup)/(m.R-(1.0-m.lambda))



    #Find the location of 0.0 in the borrowing grid
    a0Ind=searchsortedfirst(aGrid,zero(F))

    #Check that the result is actually 0.0, i.e. the grid has not been misspecified as an evenly spaced grid which skips 0.0
    @assert aGrid[a0Ind]==zero(F)


    #q=1/R*(lambda+(1-lambda)*(coup+q)
    #(R-(1-lambda))*q=lambda+(1-lambda)*coup
    #q=(lambda+(1-lambda)*coup)/(R-(1-lambda))


    netRevM0A0Grid=zeros(adim,ydim)

    #Fill the grid of consumption values when m=0 and b'=(1-lambda)*b
    for i in 1:ydim
        for j in 1:adim
            netRevM0A0Grid[j,i]=income.yGrid[i]+(m.lambda+(1.0-m.lambda)*m.coup)*aGrid[j]
        end
    end

    #Initialize the grids for consumption when m=0 and b'=(1-lambda)*b and for when m=0 and b' takes a specified value
    consM0Grid=zeros(adim,adim,ydim)

    #Fill the grid of consumption values when m=0 and b'=(1-lambda)*b
    for i in 1:(m.yParams.yPoints)
        for j in 1:(m.aPoints)
            view(consM0Grid,:,j,i).=netRevM0A0Grid[j,i].-qBase.*view(aGridIncr,:,j)
        end
    end


    #Initialize the array of maximum feasible choices of savings
    maxAPInd=ones(S,adim,ydim)

    for i in 1:ydim
        for j in 1:adim
            apInd=m.aPoints
            maxAPInd[j,i]=adim
            while apInd>1
                if consM0Grid[apInd,j,i]+income.mBounds[2]>0.0
                    maxAPInd[j,i]=apInd
                    break
                else
                    apInd+=-1
                end
            end
        end
    end





    #Generate the set of obects used to speed the grid search for the highest and lowest values of b' chosen at every (y,b) state
    apSearchParams=genBisectSearchGrids(adim)

    #construct and return the longTermBondRDSaveEval object whose components are constructed above
    outputEval=LongTermBondRDSaveEval(income,VF,pol,aGrid,aGridIncr,qBase,netRevM0A0Grid,consM0Grid,maxAPInd,a0Ind,[true],apSearchParams)


    return outputEval

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
function solveRepayChoice(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},i::S,j::S,x::F,lb::S,ub::S) where{F<:Real,S<:Integer,T}
    #If the lower and upper bound are the same, return the that index and the corresponding value as the solution
    tempC=0.0
    maxC=0.0
    if lb==ub
        tempC=s.consM0Grid[lb,j,i]+x
        return lb, u(m,tempC)+m.beta*s.VF.EVGrid[lb,i],tempC>zero(F)
    else
        #Generate temporary variables containing (in order):
        #1. the maximum consumption observed over levels of next period borrowing examined
        #2. the value of consumption for the level of next period borrowing currently being considered
        #3. the maximum available value thusfar observed for the government
        #4. the index of next period borrowing which results in the value of 3.

        #Initialize these to correspond to the values when the government borrows at the level indicated by lb
        maxC=s.consM0Grid[ub,j,i]+x
        tempC=s.consM0Grid[ub,j,i]+x
        maxVR=u(m,tempC)+m.beta*s.VF.EVGrid[ub,i]
        maxAPInd=ub

        #Iterate over the remaining values of next period borrowing
        for jj in (ub-1):-1:(lb)

            #Set the value of consumption for the level of next period borrowing currently under consideration
            tempC=s.consM0Grid[jj,j,i]+x

            #If that value is higher than the maximum value observed thusfar, calculate the government's value under this choice
            if tempC>maxC
                tempVR=u(m,tempC)+m.beta*s.VF.EVGrid[jj,i]

                #If the value of the current choice is strictly greater than the highest value thusfar observed, set the maximum value and maximizing value variables accordingly
                if tempVR>maxVR
                    maxVR=tempVR
                    maxAPInd=jj
                end

                #Track that the maximum value of consumption observed thusfar has changed
                maxC=tempC
            end
        end
        tempC=s.consM0Grid[maxAPInd,j,i]+x

        #Return the maximizing value of the next period borrowing index, the value to the government associated with it,
        #and whether the optimal choice of asset resulted in strictly positive consumption
        return maxAPInd, maxVR,tempC>zero(F)
    end
end



#solveRepayInterior! is probably the most important function in this library. It constructs the full borrowing policy functions and
#default policy functions of the government (should the government repay under any realization of m).
#Its arguments are:
#1. m: the model specification
#2. s: the list of objects used to solve the model
#3. i: the index of the persistent income state
#4. j: the index of the incoming asset state
#5. setExactP: a Boolean variable indicating whether to calculate the exact transition probabilities implied by the policy functions

function solveRepayInterior!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},i::S,j::S) where{F<:Real,S<:Integer,T}
    #Set aliases for the index of assets chosen when m takes its highest possible value and the index of assets chosen when m takes its highest possible value
    lb=s.pol.apLowM[j,i]
    ub=s.pol.apHighM[j,i]
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

        s.pol.mAPThreshold[i][1,j]=s.income.mBounds[2]
        s.pol.apPolicy[i][1,j]=lb
        #When the government never defaults in income/debt state (i,j), further set:
        #3. the first index of the borrowing policy function which is relevant for the integration to 1
        #4. the threshold for default to the lowest possible value of the m shock
        #5. the length of the list of thresholds for this state to 1

        s.pol.mListLength[j,i]=1
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
            tempCM0Old=s.consM0Grid[tempOldAP,j,i]
            tempCM0New=s.consM0Grid[tempNewAP,j,i]

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
            tempEVDiff=m.beta*(s.VF.EVGrid[tempOldAP,i]-s.VF.EVGrid[tempNewAP,i])
            if (tempCM0New>tempCM0Old)||((tempCM0New==tempCM0Old)&&(tempEVDiff<0.0))

                tempOldAP=s.pol.apPolicy[i][tempLength-1,j]
                tempLength+=-1

            else
                #Here we solve the equation 1/(C_New+m)-1/(C_Old+m)+beta*(V_Old-V_New)=0 for m. Some algebra yields:
                #(C_Old-C_New)/((C_New+m)*(C_Old+m))+beta*(V_Old-V_New)=0
                #(C_New+m)*(C_Old+m)+(C_Old-C_New)/(beta*(V_Old-V_New)=0
                #m^2+(CNew+C_Old)*m+C_New*C_Old+(C_Old-C_New)/(beta*(V_Old-V_New))=0
                #This is a parabola which opens up. We are interested in the rightmost root (the other one is, I believe, in general below whatever previous threshold was found)
                #We use the quadratic formula to do so below.

                #Calculate the term (beta*(V_Old-V_New))
                #tempEVDiff=m.beta*(s.VF.EVGrid[tempOldAP,i]-s.VF.EVGrid[tempNewAP,i])

                #tempEVDiff should, in general always be negative (if Z(y,a') is increasing in a'). If it is 0.0, set it to a negative value of very small magnitude.
                if tempEVDiff==0.0
                    tempEVDiff=-eps(F)
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
                    s.pol.mAPThreshold[i][tempLength,j]=tempMThres
                    s.pol.apPolicy[i][tempLength,j]=tempOldAP
                    tempOldAP=tempNewAP
                    tempNewAP=tempOldAP+1
                    tempLength+=1
                else
                    #If the borrowing policy function currently contains more entries than just what occurs at the lowest realization of m, we first check whether the threshold just
                    #calculated is less than the threshold above which tempOldAP SHOULD be chosen. If that is the case, then tempOldAP is NEVER chosen. We set tempOldAP to the value of
                    #the entry just above it in the borrowing policy function and deincrement the length counter (so that tempOldAP's presence in the borrowing policy function will
                    #eventually be overwritten and replaced with the proper value). Note that in this case tempNewAP is not changed.
                    if tempMThres<s.pol.mAPThreshold[i][tempLength-1,j]
                        tempOldAP=s.pol.apPolicy[i][tempLength-1,j]
                        tempLength+=-1
                    else
                        #If the new threshold is weakly greater than the old one:
                        #1. we add the upper bound of m values at which tempOldAP is chosen to the borrowing policy function
                        #2. record that it is chosen in the relevant interval
                        #3. make the new index for which an upper bound is to be determined the one which was used to find the upper bound in 1.
                        #4. set the index which we will attempt to used to find an upper bound for the m values at which the index in #3. is chosen to 1 less than the index in #3.
                        #5. increment the length counter.
                        s.pol.mAPThreshold[i][tempLength,j]=tempMThres
                        s.pol.apPolicy[i][tempLength,j]=tempOldAP
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
        s.pol.mAPThreshold[i][tempLength,j]=s.income.mBounds[2]
        s.pol.apPolicy[i][tempLength,j]=ub
        s.pol.mListLength[j,i]=tempLength
    end

    return s
end


#solveRepayRow! is one of the workhorse functions of this library. It solves the government's problem under repayment for every possible incoming level of borrowing
#given a specific value for the persistent component of income. It takes as its arguments:
#1. m: a model specification
#2. s: the collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
#4. i: the index of the value of the persistent component of income
#5. setExactP: a Boolean variable indicating whether or not to calculate exact transition probabilities

#Its output is just modified versions of s and g2

function solveRepayRow!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},g2::LongTermBondRDSaveUpdate{F},i::S) where{F<:Real,S<:Integer,T}

    tempAInd=1

    #initialize temporary variables for the upper and lower bounds of the grid search as well as the index of current asset holdings
    tempLB=1
    tempUB=s.maxAPInd[1,i]

    #Solve the government's problem at the lowest value of incoming borrowing under the lowest m shock
    s.pol.apLowM[1,i],tempVLow,g2.feasibleSolutionL[1,i]=solveRepayChoice(m,s,i,1,s.income.mBounds[1],tempLB,tempUB)

    #Mark that the problem has indeed been solved at the lowest value of incoming borrowing under the lowest m shock
    g2.solveMarkL[1,i]=true

    if g2.feasibleSolutionL[1,i]==true
        tempLB=s.pol.apLowM[1,i]
    end

    s.pol.apHighM[1,i],tempVHigh,g2.feasibleSolutionH[1,i]=solveRepayChoice(m,s,i,1,s.income.mBounds[2],tempLB,tempUB)

    #Mark that the problem has indeed been solved at the lowest value of incoming borrowing under the highest m shock
    g2.solveMarkH[1,i]=true


    #reinitialize temporary variables for the upper and lower bounds of the grid search
    tempLB=1
    tempUB=s.maxAPInd[m.aPoints,i]

    #Perform the exact same steps for the highest level of incoming borrowing
    s.pol.apLowM[m.aPoints,i],tempVLow,g2.feasibleSolutionL[m.aPoints,i]=solveRepayChoice(m,s,i,m.aPoints,s.income.mBounds[1],tempLB,tempUB)

    g2.solveMarkL[m.aPoints,i]=true


    if g2.feasibleSolutionL[m.aPoints,i]==true
        tempLB=s.pol.apLowM[m.aPoints,i]
    end

    s.pol.apHighM[m.aPoints,i],tempVHigh,g2.feasibleSolutionH[m.aPoints,i]=solveRepayChoice(m,s,i,m.aPoints,s.income.mBounds[2],tempLB,tempUB)

    g2.solveMarkH[m.aPoints,i]=true



    #Loop over the number of levels in the search parameters
    for levelInd in 1:(s.apSearchParams.numLevels)
        #Loop over the number of points in the current level
        for pointInd in 1:(s.apSearchParams.levelPoints[levelInd])
            #Set the index of the point currently under consideration
            tempAInd=s.apSearchParams.pointGrids[levelInd][pointInd]
            #If the problem has not been solved at this index for the highest realization of m, solve it. Otherwise (i.e. if the government was found to have defaulted at a
            #lower debt level under the highest realization of m), skip solving it, since the result does not matter.
            if g2.solveMarkH[tempAInd,i]==false

                #Get the indices associated with points for which the problem has already been solved whose solutions, if feasible, form upper and lower bounds
                #for the solution associated with the current index.
                tempUBInd=s.apSearchParams.pointUB[levelInd][pointInd]
                tempLBInd=s.apSearchParams.pointLB[levelInd][pointInd]

                #Check whether a feasible solution was found at those two points. If so, use the best bound available. Otherwise, use the alternative bounds already known.
                tempUB=ifelse(g2.feasibleSolutionH[tempUBInd,i]==true,min(s.maxAPInd[tempAInd,i],s.pol.apHighM[tempUBInd,i]),s.maxAPInd[tempAInd,i])
                tempLB=ifelse(g2.feasibleSolutionH[tempLBInd,i]==true,s.pol.apHighM[tempLBInd,i],1)

                #Solve the problem, check whether the solution resulted in strictly positive consumption, and mark that we have solved the problem at the current index
                s.pol.apHighM[tempAInd,i],tempVHigh,g2.feasibleSolutionH[tempAInd,i]=solveRepayChoice(m,s,i,tempAInd,s.income.mBounds[2],tempLB,tempUB)

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

            end
            #If the problem has not been solved at this index for the lowest realization of m, solve it. Otherwise (i.e. if the government was found to have defaulted at a
            #lower debt level under the highest realization of m), skip solving it, since the result does not matter.
            if g2.solveMarkL[tempAInd,i]==false
                #Get the indices associated with points for which the problem has already been solved whose solutions, if feasible, form upper and lower bounds
                #for the solution associated with the current index.
                tempUBInd=s.apSearchParams.pointUB[levelInd][pointInd]
                tempLBInd=s.apSearchParams.pointLB[levelInd][pointInd]

                #Check whether a feasible solution was found at those two points. If so, use the best bound available. Otherwise, use the alternative bounds already known.
                tempUB=ifelse(g2.feasibleSolutionL[tempUBInd,i]==true,min(s.pol.apLowM[tempUBInd,i],s.maxAPInd[tempAInd,i]),s.maxAPInd[tempAInd,i])
                tempLB=ifelse(g2.feasibleSolutionL[tempLBInd,i]==true,s.pol.apLowM[tempLBInd,i],1)

                #Solve the problem, check whether the solution resulted in strictly positive consumption, and mark that we have solved the problem at the current index
                s.pol.apLowM[tempAInd,i],tempVLow,g2.feasibleSolutionL[tempAInd,i]=solveRepayChoice(m,s,i,tempAInd,s.income.mBounds[1],tempLB,tempUB)

                g2.solveMarkL[tempAInd,i]=true

                #Check whether the value to the government under the lowest possible value of the m shock is strictly greater than the initial value of default.
                #If this is the case, then the government never defaults in this state, and we can set the neverDefault variable to true. Otherwise, there are at least some values
                #of m at which the government defaults, so we set it to false.
            end
        end
    end

    #Loop over levels of debt.
    for j in 1:m.aPoints

        solveRepayInterior!(m,s,i,j)

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
function integrateMVAlmostExact(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},i::S,j::S) where{F<:Real,S<:Integer,T}
    outValV=0.0

    #If the government never defaults, then we begin the integration over the result under repayment at the first interval.
    #We set the location of the current upper bound in the m grid to the second point, the index of the relevant borrowing policy entry to the first,
    #and the most recent value of m to its lower bound
    if i!=m.yParams.yPoints
        mUBInd=2
        aPolicyInd=1
        lastMVal=s.income.mGrid[1]
        tempCBase=0.0
        tempCHigh=0.0
        tempCLow=0.0

        #Loop until we exit the upper boundary of the m space (the second check is to ensure that we never exit the range of the current borrowing policy function;
        #it should never matter, and if it does, the process should in general not lead to convergence).
        while (mUBInd<=(m.mParams.mPoints))&&(aPolicyInd<=(s.pol.mListLength[j,i]))
            #If the upper bound for the range in which the current entry of the borrowing policy function is valid is strictly less than the upper bound of the current
            #m interval, add the relevant contribution of that entry, update the last value of m reached, and increment only the variable marking our position in the borrowing
            #policy function
            if s.pol.mAPThreshold[i][aPolicyInd,j]<(s.income.mGrid[mUBInd])
                tempCBase=s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]
                tempCHigh=tempCBase+s.pol.mAPThreshold[i][aPolicyInd,j]
                tempCLow=tempCBase+lastMVal
                outValV+=s.income.mProb[mUBInd-1]*((log(tempCLow)-log(tempCHigh))/s.income.mRes+(s.pol.mAPThreshold[i][aPolicyInd,j]-lastMVal)/s.income.mRes*m.beta*s.VF.EVGrid[s.pol.apPolicy[i][aPolicyInd,j],i])

                lastMVal=s.pol.mAPThreshold[i][aPolicyInd,j]
                aPolicyInd+=1
            else
                #Otherwise, add the relevant contribution of the current entry which lies in this m interval, update the last value of m reached to be its upper bound,
                #and then increment only the variable marking our position in the m grid
                tempCBase=s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i]
                tempCHigh=tempCBase+s.income.mGrid[mUBInd]
                tempCLow=tempCBase+lastMVal
                outValV+=s.income.mProb[mUBInd-1]*((log(tempCLow)-log(tempCHigh))/s.income.mRes+(s.income.mGrid[mUBInd]-lastMVal)/s.income.mRes*m.beta*s.VF.EVGrid[s.pol.apPolicy[i][aPolicyInd,j],i])

                lastMVal=s.income.mGrid[mUBInd]
                mUBInd+=1
            end
        end

    else
        aPolicyInd=searchsortedfirst(s.pol.mAPThreshold[i][1:s.pol.mListLength[j,i],j],0.0)
        outValV+=u(m,s.consM0Grid[s.pol.apPolicy[i][aPolicyInd,j],j,i])+m.beta*s.VF.EVGrid[s.pol.apPolicy[i][aPolicyInd,j],i]
    end

    return outValV
end



#updateV! performs the full integration across states to calculate V(y,b)=E_m[W(y,m,b)]. Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
#4. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of vNew
function updateV!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},g2::LongTermBondRDSaveUpdate{F})  where{F<:Real,S<:Integer,T}

    Threads.@threads for i in 1:(m.yParams.yPoints)
    #for i in 1:(m.yParams.yPoints)
        for j in 1:(m.aPoints)
            g2.VF.vGrid[j,i]=integrateMVAlmostExact(m,s,i,j)
        end
    end


    return g2
end

#updateEV! performs the second integration step in updating the main value function. Specifically, it calculates Z(y,b')=E_y'[V(y',b')].
#Its arguments are:
#1. m: a model specification
#2. s: a collection of objects used to solve the model
#3. g2: the collection of objects used, along with the policy functions in s, to solve the model. update the guesses of q and Z
function updateEV!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},g2::LongTermBondRDSaveUpdate{F})  where{F<:Real,S<:Integer,T}
    mul!(g2.VF.EVGrid,g2.VF.vGrid,s.income.yTMat)
    return g2
end


#vfiGOneStep! is the main function of the for solving method for the model in Chatterjee and Eyigungor (2012)
#1. m: a model specification
#2. s: a collection of objects used to solve the model which will be modified
#3. tol: the tolerance for convergence for both the value functions and the bond price function
#4. maxIter: the maximum number of iterations
#5. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of s
function vfiGOneStep!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},tol::F,maxIter::S) where{F<:Real,S<:Integer,T}




    #Initialize new copies of the value functions and the bond price function
    g2=makeLTBUpdate(m,s)


    #Initialize  a counter for the number of iterations and measures of sup norm distance between various objects
    iCount=1
    vFDist=tol+1.0
    EVDist=tol+1.0

    vDist=tol+1.0
    maxDist=tol+1.0

    #Iterate until the maximum number of iterations has been reached or the sup norm distance between successive iterations drops below the tolerance for convergence
    while (iCount<=maxIter)&(maxDist>=tol)

        #Reset all the objects used to speed up the grid search step
        g2.solveMarkL.=false
        g2.solveMarkH.=false
        g2.feasibleSolutionL.=false
        g2.feasibleSolutionH.=false


        #Iterate over states of the persistent component of income
        Threads.@threads for i in 1:(m.yParams.yPoints)
        #for i in 1:(m.yParams.yPoints)
            #solve the problem at the current state of the persistent component of income
            solveRepayRow!(m,s,g2,i)
        end

        #Perform the integration steps to obtain new versions of V, Z, and q
        updateV!(m,s,g2)
        updateEV!(m,s,g2)


        #Calculate the measures of distance between successive iterations for V, Z, and q
        vFDist=maximum(abs.(g2.VF.vGrid-s.VF.vGrid))
        EVDist=maximum(abs.(g2.VF.EVGrid-s.VF.EVGrid))

        #Update the guesses for V, Z, and q
        s.VF.vGrid.=g2.VF.vGrid
        s.VF.EVGrid.=g2.VF.EVGrid



        #Calculate the maximum of the various distance measures for the value functions
        vDist=max(vFDist,EVDist)

        #Calculate the maximum of all the various distance measures
        maxDist=vDist

        #At each iteration which is a multiple of 5% of the maximum number of iterations, print the current iteration number and variables tracking the convergence criteria
        if mod(iCount,100)==0
            # println([iCount,EVDist,vFDist])
            print(".")
        end

        #Increment the iteration counter
        iCount+=1


    end



    #Print final set of convergence statistics
    # println("\n", [iCount-1,EVDist,vFDist])

    #Update the marker tracking convergence of the value functions
    if vDist>tol
        s.noConvergeMark[1]=true
    else
        s.noConvergeMark[1]=false
    end


    #Return the modified collection of objects used to solve the model
    return s
end


#vfiGOneStep! is the main function of the for solving method for the model in Chatterjee and Eyigungor (2012)
#1. m: a model specification
#2. s: a collection of objects used to solve the model which will be modified
#3. tol: the tolerance for convergence for both the value functions and the bond price function
#4. maxIter: the maximum number of iterations
#5. vIntExact: a true/false variable indicating whether quadrature or midpoint integration is to be used. True corresponds to quadrature

#Its output is simply a modified version of s
function vfiGOneStep!(m::LongTermBondRDSaveSpec{F,S},s::LongTermBondRDSaveEval{F,S,T},tol::F,maxIter::S,mixFacV::F) where{F<:Real,S<:Integer,T}




    #Initialize new copies of the value functions and the bond price function
    g2=makeLTBUpdate(m,s)




    #Initialize  a counter for the number of iterations and measures of sup norm distance between various objects
    iCount=1
    vFDist=tol+1.0
    EVDist=tol+1.0

    vDist=tol+1.0
    maxDist=tol+1.0

    #Iterate until the maximum number of iterations has been reached or the sup norm distance between successive iterations drops below the tolerance for convergence
    while (iCount<=maxIter)&(maxDist>=tol)

        #Reset all the objects used to speed up the grid search step
        g2.solveMarkL.=false
        g2.solveMarkH.=false
        g2.feasibleSolutionL.=false
        g2.feasibleSolutionH.=false


        #Iterate over states of the persistent component of income
        Threads.@threads for i in 1:(m.yParams.yPoints)
        #for i in 1:(m.yParams.yPoints)
            #solve the problem at the current state of the persistent component of income
            solveRepayRow!(m,s,g2,i)
        end

        #Perform the integration steps to obtain new versions of V, Z, and q
        updateV!(m,s,g2)
        updateEV!(m,s,g2)


        #Calculate the measures of distance between successive iterations for V, Z, and q
        vFDist=maximum(abs.(g2.VF.vGrid-s.VF.vGrid))
        EVDist=maximum(abs.(g2.VF.EVGrid-s.VF.EVGrid))

        #Update the guesses for V, Z, and q
        s.VF.vGrid.=mixFacV*s.VF.vGrid.+(1.0-mixFacV)*g2.VF.vGrid
        s.VF.EVGrid.=mixFacV*s.VF.EVGrid.+(1.0-mixFacV)*g2.VF.EVGrid



        #Calculate the maximum of the various distance measures for the value functions
        vDist=max(vFDist,EVDist)

        #Calculate the maximum of all the various distance measures
        maxDist=vDist

        #At each iteration which is a multiple of 5% of the maximum number of iterations, print the current iteration number and variables tracking the convergence criteria
        if mod(iCount,100)==0
            # println([iCount,EVDist,vFDist])
            print(".")
        end

        #Increment the iteration counter
        iCount+=1


    end



    #Print final set of convergence statistics
    # println("\n", [iCount-1,EVDist,vFDist])

    #Update the marker tracking convergence of the value functions
    if vDist>tol
        s.noConvergeMark[1]=true
    else
        s.noConvergeMark[1]=false
    end


    #Return the modified collection of objects used to solve the model
    return s
end


function writecsv(fileName::String,gridToWrite::Array{F,N}) where{F,N}
    writedlm(fileName,gridToWrite,',')
end


function readcsv(fileName::String)
    readdlm(fileName,',')
end

#writeEssentials is a Quick function for saving the pieces of a model necessary for producing figures as well as the full set of objects
#which can be used to reconstruct the computed model given the specification only.
function writeEssentials(s::LongTermBondRDSaveEval{F,S,T},filePrefix::String,fileDir::String) where{F<:Real,S<:Integer,T}
    writecsv(
        joinpath(fileDir, filePrefix*"_yGrid.csv"),
        s.income.yGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_aGrid.csv"),
        s.aGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_vGrid.csv"),
        s.VF.vGrid
    )
    writecsv(
        joinpath(fileDir, filePrefix*"_EVGrid.csv"),
        s.VF.EVGrid
    )
end
