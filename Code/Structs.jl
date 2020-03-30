#######################################################################
# Income structs 
#######################################################################


#ar1RDParams is a parametrically typed immutable containing:
#1. yPoints: the number of points to be used in discretizing the AR(1) Process
#2. rho: the persistence parameter of the process
#3. eta2: the variance of the innovations of the process
#4. mu: the permanent component of the process
#5. stdSpan: the number of standard deviations to each side of the mean which should fall within the grid
#6. inflateEndpoints: a Boolean variable marking whether the probability of jumping to a point outside the grid should be assigned to the closest endpoint (true)
#or discarded (effectively reassigning it to the remaining points in proportion to their probability of being reached)
#Throughout, we use the convention y_t=mu+rho*y_t-1+sqrt(eta2)*e_t

struct AR1RDParams{F<:Real,S<:Integer}
    yPoints::S
    rho::F
    eta2::F
    muVec::Array{F,1}
    muTMat::Array{F,2}
    stdSpan::Array{F,1}
    inflateEndpoints::Bool
    ergodicReentry::Bool
end


#iidParams is a parametrically typed immutable containing:
#1. mPoints: the number of points to be used in discretizing the AR(1) Process
#2. eta2: the variance of the process
#3. mu: the mean of the process
#4. stdSpan: the number of standard deviations to each side of the mean which should fall within the grid
struct IIDParams{F<:Real,S<:Integer}
    mPoints::S
    epsilon2::F
    mu::F
    stdSpan::F
end

#######################################################################
# Grid and Distribution structs 
#######################################################################


#gridSearchParams is a parametrically typed immutable containing the index information required to quickly access point and bound indices to perform a bisection style grid search which exploits policy function monotonicity
#Its contents are:
#1. boundGrids: a one dimensional array containing the one dimensional arrays of bound indices at each level
#2. pointGrids:a one dimensional array containing the one dimensional arrays of point indices at each level
#3. pointLB: a one dimensional array containing the one dimensional arrays of lower bound indices for each point at each level
#4. pointUB: a one dimensional array containing the one dimensional arrays of upper bound indices for each point at each level
#5. boundGridPoints: a one dimensional array containing the length of the arrays in boundGrids
#6. levelPoints: a one dimensional array containing the length of the arrays in pointGrids (and therefore also in pointLB and pointUB)
#7. numLevels: the length all the above arrays (at the top level, i.e. the number of actual arrays contained in boundGrids rather than the number of points in all the arrays it contains)
struct GridSearchParams{S<:Integer}
    boundGrids::Array{Array{S,1},1}
    pointGrids::Array{Array{S,1},1}
    pointLB::Array{Array{S,1},1}
    pointUB::Array{Array{S,1},1}
    boundGridPoints::Array{S,1}
    levelPoints::Array{S,1}
    numLevels::S
end


#debtDist is a parametrically typed immutable describing the stationary, ergodic joint distirbution of default state, persistent income state, and current borrowing state. Its contents are:

#1. repayDist: a two dimensional array (axes are 1. incoming borrowing, 2. persistent income state) containing the noninflated joint distribution of income and borrowing, conditional
#on the government having access to financial markets (where by noninflated, we mean that we have not divided the conditional distribution by theprobability of the condition being true)
#2. defaultDist: a one dimensional array (axes are 1. income) containing the noninflated distribution of income, conditional on the government being in default
struct DebtDist{F<:Real}
    repayDist::Array{F,2}
    defaultDist::Array{F,1}
end


#######################################################################
# Borrowing equilibrium structs 
#######################################################################


#longTermBondRDBorrowSpec is a parametrically typed immutable containing enough of the information required to specify the model in Chatterjee and Eyigungor (2012)
#that the remainder can be fully deduced and constructed. Its contents are:

#1. beta: the government's discount factor
#2. theta: the probability of reentry to international markets when in financial autarky
#3. gamma: the coefficient of relative risk aversion for CRRA preferences (which we assume both the government and consumer have, throughout)
#4. hpen0: either the h in y_def(s)=min(y(s),h*E[y(s)]), i.e. the output penalty when in default, or d0 in y_def(s)=y(s)-max(0,d0*y(s)+d1*y(s)^2)
#5. hpen1: d1 in y_def(s)=y(s)-max(0,d0*y(s)+d1*y(s)^2)
#6. R: the gross interational interest rate
#7. lambda: the maturity parameter of the debt structure
#8. coup: the coupon parameter of the debt structure
#9. aBounds: the boundaries of the borrowing grid [min, max]
#10. aPoints: the number of points to be used in the borrowing grid
#11. yParams: the specification of the persistent component of the income process and how it shall be approximated
#12. iidParams: the specification of the iid component of the income process
#13. simplePen: a boolean variable which is true if y_def(s)=min(y(s),h*E[y(s)]) and false if y_def(s)=y(s)-max(0,d0*y(s)+d1*y(s)^2)
#14. mixFacQ: the convex combination parameter in q^k+1 = mixFacQ*q^k + (1-mixFacQ)*H(q^k) where H is the operator which given government policiesBorrow conditional on q^k, updates the bond price function
#15. mixFacV: the convex combination parameter in Z^k+1 = mixFacV*Z^k + (1-mixFacQ)*T(Z^k) where T is the operator which given government policiesBorrow conditional on q^k, updates the expected continuation value function

#Throughout the code, the variable m is always an object of type longTermBondRDBorrowSpec

struct LongTermBondRDBorrowSpec{F<:Real,S<:Integer}
    beta::F
    theta::F
    gamma::F
    R::F
    lambda::F
    coup::F
    hpen0::F
    hpen1::F
    aBounds::Array{F,1}
    aPoints::S
    yParams::AR1RDParams{F,S}
    mParams::IIDParams{F,S}
    simplePen::Bool
    mixFacQ::F
    mixFacV::F
end

#vFuncBorrow is a parametrically typed immutable containing the flow utilities in the period of default and thereafter and
#the value function/continuation value function grids. Its contents are:
#1. vDFutFlow: a one dimensional array (axes are 1. current persistent income) containing the the expected value of the
#flow utility term conditional on the realization of the persistent component while in default
#2. vDInitFlow: a one dimensional array (axes are 1. current persistent income) containing the the value of the
#flow utility term conditional on the realization of the persistent component when the government enters the period
#with access to financial markets and then defaults
#3. vGrid: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing
#the expected values of the value function conditional on the realization of y and a, i.e. V(y,a)=E_m[W(y,m,a)]
#4. EVGrid: a two dimensional array (axes are 1. next period asset level and 2. current persistent income) containing
#the grid of the expected continuation value function Z(y,a')=E[V(y',a')|y]
#5. vDInitGrid: a one dimensional array (axes are 1. current persistent income) containing the value to the government of defaulting
#when it enters the period with access to financial markets.
#6. vDFutGrid: a one dimensional array (axes are 1. current persistent income) containing the value to the government
#when it enters the period without access to financial markets (i.e. already in default)
#7. EVDGrid: a one dimensional array (axes are 1. current persistent income) containing
#the grid of the expected continuation value function Z^D(y)=E[theta*V(y',0)+(1.0-theta)*V^D(y')|y]

struct VFuncBorrow{F<:Real}
    vDFutFlow::Array{F,1}
    vDInitFlow::Array{F,1}
    vDFutGrid::Array{F,1}
    vDInitGrid::Array{F,1}
    EVDGrid::Array{F,1}
    vGrid::Array{F,2}
    EVGrid::Array{F,2}
end

#policiesBorrow is a parametrically typed immutable containing the government's policy functions. Its contents are:
#1. apLowM: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the government's borrowing policy, conditional on repayment,
#when the m shock takes its lowest value.
#2. apHighM: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the government's borrowing policy, conditional on repayment,
#when the m shock takes its highest value.

#3., 4., 5., and 6.  are mListLength, mAPThreshold, apPolicy, apProbability
#mListLength, mAPThreshold, and apPolicy fully describe the government's policy function conditional on repayment.

#7., 8., and 9. are defThreshold, defProb, and firstRepayInd
#When the above mentioned three three are combined additionally with defThreshold and firstRepayInd, they completely describe the government's policies,
#except in two special cases described just below

#10. and 11. are alwaysDefault and neverDefault
#these arrays contain marker variables for whether we can take shortcuts during the steps involving integration over m.

#apProbability and defProb are included for convenience. Combined with apPolicy and mListLength, they fully describe the
#exact transition probabilities implied by government policies.

#Specifically these objects are:

#3. mListLength: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the number of distinct
#levels of asset which, conditional on repayment, the government might choose depending on the value of m. Since the rest of the government's
#policies conditional on repayment is stored in sparse matrices which are NOT zeroed at irrelevant values once the government's problem is
#solved again, we have to keep track of exactly how many cells are actually in use

#4. mAPThreshold: a one dimensional array (axes are 1. current persistent income) containing two dimensional sparse matrices.
#Each of these matrices is square NxN with N equal to the number of borrowing points (so that each column can potentially hold the longest
#possible policy function). For each of these matrices, the second axis is current asset level. Given a specific matrix and column, each cell
#from row 1 to the row indicated by the corresponding entry of mListLength indicates the upper bound of the values of m for which the government,
#conditional on repayment, chooses the borrowing policy specified in apPolicy.

#5. apPolicy: a one dimensional array (axes are 1. current persistent income) containing two dimensional sparse matrices.
#The layout of each is exactly as for those in mAPThreshold. Given the matrix for persistent income level i and column for asset level j,
#the entries of apPolicy[i][:,j] in rows 1 to the value of mListLength[j,i] indicate the policy choices associated with the thresholds in mAPThreshold.

#Combining these three allows us to say that, when persistent income level i and current asset level j, conditional on repayment,
#for m between it minimum value and the mAPThreshold[i][1,j], the government chooses apPolicy[i][1,j], and if mListLength[j,i] is strictly greater than 1,
#then it chooses apPolicy[i][jj,j] for m in between mAPThreshold[i][jj-1,j] and mAPThreshold[i][jj,j] for jj between 2 and mListLength[j,i], inclusive.

#6. apProbability: a one dimensional array (axes are 1. current persistent income) containing two dimensional sparse matrices.
#The layout of each is exactly as for those in mAPThreshold and apPolicy. It contains the probability that apPolicy[i][jj,j] is chosen
#for jj between 1 and mListLength[j,i], inclusive. These probabilities are NOT conditional on repayment. They are calculated taking into account
#the threshold level of m below which default occurs.

#7. defThreshold: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the levels of m below
#which default occurs.

#8. defProbability: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the probabilities of default.

#9. firstRepayInd: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the first indices of mAPThreshold,
#apPolicy, and apProbability which is relevant for integrating the government's value function over realizations of m.

#10. alwaysDefault: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing Boolean variables indicating
#whether the government always defaults at each persistent income and asset state.

#11. neverDefault: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing Boolean variables indicating
#whether the government never defaults at each persistent income and asset state.

struct PoliciesBorrow{F<:Real,S<:Integer}
    apLowM::Array{S,2}
    apHighM::Array{S,2}
    mListLength::Array{S,2}
    mAPThreshold::Array{SparseMatrixCSC{F,S},1}
    apPolicy::Array{SparseMatrixCSC{S,S},1}
    firstRepayInd::Array{S,2}
    defThreshold::Array{F,2}
    alwaysDefault::Array{Bool,2}
    neverDefault::Array{Bool,2}
end



#combinedIncomeProcessBorrow is a parametrically typed immutable which contains the full specification of the approximation (if any)
#of both components of the income process. Its contents are:
#1. yGrid: a one dimensional array containing the persistent component of the income process when the country is not in financial autarky
#2. yDefGrid: a one dimemsional array containing the persistent component of the income process when the country is in financial autarky
#3. yTMat: the transition matrix for the persistent component of the income process. It is produced using the convention that columns sum to 1.0
#4. yMean: the long term mean value (analytically) of the peristent component of the income process.
#5. mBounds: the minimum and maximum values assumed by the iid component of the income process.
#6. mGrid: a one dimensional array containing the grid of points which define the intervals for the m shock to be used in the integration steps
#7. mMidPoints: a one dimensional array containing the midpoints of the intervals in 6.
#8. mProb: the probability of realizing m in each interval defined in 6.
#9. mRes: the resolution of the grid in 6.
#10. mInfl: the inflation factor for m shock probabilities, i.e. 1/(CDF(m_max)-CDF(m_min))
#11. mDist: the distribution of the m shock (will be an object generated by the Distributions package)
struct CombinedIncomeProcessBorrow{F<:Real,T}
    yGrid::Array{F,1}
    yDefGrid::Array{F,1}
    yTMat::Array{F,2}
    yMean::F
    mBounds::Array{F,1}
    mGrid::Array{F,1}
    mMidPoints::Array{F,1}
    mProb::Array{F,1}
    mRes::F
    mInfl::F
    mDist::T
end

#longTermBondRDBorrowEval is a parametrically typed immutable which provides, in combination with a longTermBondRDBorrowSpec, all the objects and information required to solve
#the model of Chatterjee and Eyigungor (2008). Its contents are:
#1. income: a combined income process object described above.
#2. VF: a value function object described above.
#3. pol: a policy function object described above.
#4. aGrid: a one dimensional array containing the values of asset which the government can choose.
#5. qGrid: a two dimensional array (axes are 1. next period asset level and 2. current persistent income) containing the bond price function
#6. netRevM0A0Grid: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the level of consumption
#when the value of m is exactly zero and the government does not raise any new debt or buy back any old debt
#7. consM0Grid: a two dimensional array (axes are 1. next period asset level, 2. current asset level and 3. current persistent income level) containing the level of consumption
#when the value of m is exactly zero and the government sets a specified new level of (axes are 1. current asset level and 2. current persistent income)
#8. maxAPInd: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the minimum index in the next period
#borrowing grid which results in strictly positive consumption under the maximum realization of m
#9. a0Ind: the index of 0.0 in the grid of (axes are 1. current asset level and 2. current persistent income) levels.
#10. noConvergeMark: a one dimensional array which simply tracks whether the model converged to the specified tolerance
#11. apSearchParams: the set of objects required for exploiting the binary monotonicity of the borrowing policy functions when finding the
#minimum and maximum values it takes for each level of income and assets

struct LongTermBondRDBorrowEval{F<:Real,S<:Integer,T}
    income::CombinedIncomeProcessBorrow{F,T}
    VF::VFuncBorrow{F}
    pol::PoliciesBorrow{F,S}
    aGrid::Array{F,1}
    aGridIncr::Array{F,2}
    qGrid::Array{F,2}
    netRevM0A0Grid::Array{F,2}
    consM0Grid::Array{F,3}
    maxAPInd::Array{S,2}
    a0Ind::S
    noConvergeMark::Array{Bool,1}
    apSearchParams::GridSearchParams{S}
end


#longTermBondRDBorrowUpdate is a parametrically typed immutable which contains the objects necessary to help solve
#the model and then calculate updated versions of the bond price function and the continuation value functions.
#Its contents are:
#1. VF: a value function object described above
#2. EVA0: a one dimensional array containing the reentry continuation value
#3. qGrid: a two dimensional array (axes are 1. next period asset level and 2. current persistent income) containing the updated bond price function
#4. qSum: a two dimensional array (axes are 1. next period asset level and 2. nexy period persistent income) containing the the value of the RHS terms
#in the functional equation describing the bond price function, i.e. qGrid[i,j]=E[1/R*qSum[i',j]|i]
#5. solveMarkH: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the highest level of the m shock has been solved during the current iteration
#6. solveMarkL: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the lowest level of the m shock has been solved during the current iteration.
#7. feasibleSolutionH: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the highest level of the m shock has been solved during the current iteration and that solution resulted in strictly positive consumption.
#8. feasibleSolutionL: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the lowest level of the m shock has been solved during the current iteration and that solution resulted in strictly positive consumption.
#9. maxAlwaysDefInd: a one dimensional array (axes are 1. current persistent income) containing the maximum asset level
#at which the government defaults under all possible realizations of the m shock

struct LongTermBondRDBorrowUpdate{F<:Real,S<:Integer}
    VF::VFuncBorrow{F}
    EVA0::Array{F,1}
    qGrid::Array{F,2}
    qSum::Array{F,2}
    solveMarkH::Array{Bool,2}
    solveMarkL::Array{Bool,2}
    feasibleSolutionH::Array{Bool,2}
    feasibleSolutionL::Array{Bool,2}
    maxAlwaysDefInd::Array{S,1}
end


#######################################################################
# Savings equilibrium structs 
#######################################################################

#longTermBondRDSaveSpec is a parametrically typed immutable containing enough of the information required to specify the model in Chatterjee and Eyigungor (2012)
#that the remainder can be fully deduced and constructed. Its contents are:

#1. beta: the government's discount factor
#2. theta: the probability of reentry to international markets when in financial autarky
#3. gamma: the coefficient of relative risk aversion for CRRA preferences (which we assume both the government and consumer have, throughout)
#4. hpen0: either the h in y_def(s)=min(y(s),h*E[y(s)]), i.e. the output penalty when in default, or d0 in y_def(s)=y(s)-max(0,d0*y(s)+d1*y(s)^2)
#5. hpen1: d1 in y_def(s)=y(s)-max(0,d0*y(s)+d1*y(s)^2)
#6. R: the gross interational interest rate
#7. lambda: the maturity parameter of the debt structure
#8. coup: the coupon parameter of the debt structure
#9. aBounds: the boundaries of the borrowing grid [min, max]
#10. aPoints: the number of points to be used in the borrowing grid
#11. yParams: the specification of the persistent component of the income process and how it shall be approximated
#12. iidParams: the specification of the iid component of the income process
#13. simplePen: a boolean variable which is true if y_def(s)=min(y(s),h*E[y(s)]) and false if y_def(s)=y(s)-max(0,d0*y(s)+d1*y(s)^2)
#14. mixFacQ: the convex combination parameter in q^k+1 = mixFacQ*q^k + (1-mixFacQ)*H(q^k) where H is the operator which given government policies conditional on q^k, updates the bond price function
#15. mixFacV: the convex combination parameter in Z^k+1 = mixFacV*Z^k + (1-mixFacQ)*T(Z^k) where T is the operator which given government policies conditional on q^k, updates the expected continuation value function

#Throughout the code, the variable m is always an object of type longTermBondRDSaveSpec

struct LongTermBondRDSaveSpec{F<:Real,S<:Integer}
    beta::F
    gamma::F
    R::F
    lambda::F
    coup::F
    aBounds::Array{F,1}
    aPoints::S
    yParams::AR1RDParams{F,S}
    mParams::IIDParams{F,S}
end

#vFunc is a parametrically typed immutable containing the flow utilities in the period of default and thereafter and
#the value function/continuation value function grids. Its contents are:
#1. vDFutFlow: a one dimensional array (axes are 1. current persistent income) containing the the expected value of the
#flow utility term conditional on the realization of the persistent component while in default
#2. vDInitFlow: a one dimensional array (axes are 1. current persistent income) containing the the value of the
#flow utility term conditional on the realization of the persistent component when the government enters the period
#with access to financial markets and then defaults
#3. vGrid: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing
#the expected values of the value function conditional on the realization of y and a, i.e. V(y,a)=E_m[W(y,m,a)]
#4. EVGrid: a two dimensional array (axes are 1. next period asset level and 2. current persistent income) containing
#the grid of the expected continuation value function Z(y,a')=E[V(y',a')|y]
#5. vDInitGrid: a one dimensional array (axes are 1. current persistent income) containing the value to the government of defaulting
#when it enters the period with access to financial markets.
#6. vDFutGrid: a one dimensional array (axes are 1. current persistent income) containing the value to the government
#when it enters the period without access to financial markets (i.e. already in default)
#7. EVDGrid: a one dimensional array (axes are 1. current persistent income) containing
#the grid of the expected continuation value function Z^D(y)=E[theta*V(y',0)+(1.0-theta)*V^D(y')|y]

struct VFuncSave{F<:Real}
    vGrid::Array{F,2}
    EVGrid::Array{F,2}
end
#policies is a parametrically typed immutable containing the government's policy functions. Its contents are:
#1. apLowM: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the government's borrowing policy, conditional on repayment,
#when the m shock takes its lowest value.
#2. apHighM: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the government's borrowing policy, conditional on repayment,
#when the m shock takes its highest value.

#3., 4., 5., and 6.  are mListLength, mAPThreshold, apPolicy, apProbability
#mListLength, mAPThreshold, and apPolicy fully describe the government's policy function conditional on repayment.

#7., 8., and 9. are defThreshold, defProb, and firstRepayInd
#When the above mentioned three three are combined additionally with defThreshold and firstRepayInd, they completely describe the government's policies,
#except in two special cases described just below

#10. and 11. are alwaysDefault and neverDefault
#these arrays contain marker variables for whether we can take shortcuts during the steps involving integration over m.

#apProbability and defProb are included for convenience. Combined with apPolicy and mListLength, they fully describe the
#exact transition probabilities implied by government policies.

#Specifically these objects are:

#3. mListLength: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the number of distinct
#levels of asset which, conditional on repayment, the government might choose depending on the value of m. Since the rest of the government's
#policies conditional on repayment is stored in sparse matrices which are NOT zeroed at irrelevant values once the government's problem is
#solved again, we have to keep track of exactly how many cells are actually in use

#4. mAPThreshold: a one dimensional array (axes are 1. current persistent income) containing two dimensional sparse matrices.
#Each of these matrices is square NxN with N equal to the number of borrowing points (so that each column can potentially hold the longest
#possible policy function). For each of these matrices, the second axis is current asset level. Given a specific matrix and column, each cell
#from row 1 to the row indicated by the corresponding entry of mListLength indicates the upper bound of the values of m for which the government,
#conditional on repayment, chooses the borrowing policy specified in apPolicy.

#5. apPolicy: a one dimensional array (axes are 1. current persistent income) containing two dimensional sparse matrices.
#The layout of each is exactly as for those in mAPThreshold. Given the matrix for persistent income level i and column for asset level j,
#the entries of apPolicy[i][:,j] in rows 1 to the value of mListLength[j,i] indicate the policy choices associated with the thresholds in mAPThreshold.

#Combining these three allows us to say that, when persistent income level i and current asset level j, conditional on repayment,
#for m between it minimum value and the mAPThreshold[i][1,j], the government chooses apPolicy[i][1,j], and if mListLength[j,i] is strictly greater than 1,
#then it chooses apPolicy[i][jj,j] for m in between mAPThreshold[i][jj-1,j] and mAPThreshold[i][jj,j] for jj between 2 and mListLength[j,i], inclusive.

#6. apProbability: a one dimensional array (axes are 1. current persistent income) containing two dimensional sparse matrices.
#The layout of each is exactly as for those in mAPThreshold and apPolicy. It contains the probability that apPolicy[i][jj,j] is chosen
#for jj between 1 and mListLength[j,i], inclusive. These probabilities are NOT conditional on repayment. They are calculated taking into account
#the threshold level of m below which default occurs.

#7. defThreshold: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the levels of m below
#which default occurs.

#8. defProbability: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the probabilities of default.

#9. firstRepayInd: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the first indices of mAPThreshold,
#apPolicy, and apProbability which is relevant for integrating the government's value function over realizations of m.

#10. alwaysDefault: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing Boolean variables indicating
#whether the government always defaults at each persistent income and asset state.

#11. neverDefault: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing Boolean variables indicating
#whether the government never defaults at each persistent income and asset state.

struct PoliciesSave{F<:Real,S<:Integer}
    apLowM::Array{S,2}
    apHighM::Array{S,2}
    mListLength::Array{S,2}
    mAPThreshold::Array{SparseMatrixCSC{F,S},1}
    apPolicy::Array{SparseMatrixCSC{S,S},1}
end

#combinedIncomeProcessSave is a parametrically typed immutable which contains the full specification of the approximation (if any)
#of both components of the income process. Its contents are:
#1. yGrid: a one dimensional array containing the persistent component of the income process when the country is not in financial autarky
#2. yDefGrid: a one dimemsional array containing the persistent component of the income process when the country is in financial autarky
#3. yTMat: the transition matrix for the persistent component of the income process. It is produced using the convention that columns sum to 1.0
#4. yMean: the long term mean value (analytically) of the peristent component of the income process.
#5. mBounds: the minimum and maximum values assumed by the iid component of the income process.
#6. mGrid: a one dimensional array containing the grid of points which define the intervals for the m shock to be used in the integration steps
#7. mMidPoints: a one dimensional array containing the midpoints of the intervals in 6.
#8. mProb: the probability of realizing m in each interval defined in 6.
#9. mRes: the resolution of the grid in 6.
#10. mInfl: the inflation factor for m shock probabilities, i.e. 1/(CDF(m_max)-CDF(m_min))
#11. mDist: the distribution of the m shock (will be an object generated by the Distributions package)
struct CombinedIncomeProcessSave{F<:Real,T}
    yGrid::Array{F,1}
    yTMat::Array{F,2}
    yMean::F
    mBounds::Array{F,1}
    mGrid::Array{F,1}
    mMidPoints::Array{F,1}
    mProb::Array{F,1}
    mRes::F
    mInfl::F
    mDist::T
end

#longTermBondRDSaveEval is a parametrically typed immutable which provides, in combination with a longTermBondRDSaveSpec, all the objects and information required to solve
#the model of Chatterjee and Eyigungor (2008). Its contents are:
#1. income: a combined income process object described above.
#2. VF: a value function object described above.
#3. pol: a policy function object described above.
#4. aGrid: a one dimensional array containing the values of asset which the government can choose.
#5. qGrid: a two dimensional array (axes are 1. next period asset level and 2. current persistent income) containing the bond price function
#6. netRevM0A0Grid: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the level of consumption
#when the value of m is exactly zero and the government does not raise any new debt or buy back any old debt
#7. consM0Grid: a two dimensional array (axes are 1. next period asset level, 2. current asset level and 3. current persistent income level) containing the level of consumption
#when the value of m is exactly zero and the government sets a specified new level of (axes are 1. current asset level and 2. current persistent income)
#8. maxAPInd: a two dimensional array (axes are 1. current asset level and 2. current persistent income) containing the minimum index in the next period
#borrowing grid which results in strictly positive consumption under the maximum realization of m
#9. a0Ind: the index of 0.0 in the grid of (axes are 1. current asset level and 2. current persistent income) levels.
#10. noConvergeMark: a one dimensional array which simply tracks whether the model converged to the specified tolerance
#11. apSearchParams: the set of objects required for exploiting the binary monotonicity of the borrowing policy functions when finding the
#minimum and maximum values it takes for each level of income and assets

struct LongTermBondRDSaveEval{F<:Real,S<:Integer,T}
    income::CombinedIncomeProcessSave{F,T}
    VF::VFuncSave{F}
    pol::PoliciesSave{F,S}
    aGrid::Array{F,1}
    aGridIncr::Array{F,2}
    qBase::F
    netRevM0A0Grid::Array{F,2}
    consM0Grid::Array{F,3}
    maxAPInd::Array{S,2}
    a0Ind::S
    noConvergeMark::Array{Bool,1}
    apSearchParams::GridSearchParams{S}
end


#longTermBondRDSaveUpdate is a parametrically typed immutable which contains the objects necessary to help solve
#the model and then calculate updated versions of the bond price function and the continuation value functions.
#Its contents are:
#1. VF: a value function object described above
#2. EVA0: a one dimensional array containing the reentry continuation value
#3. qGrid: a two dimensional array (axes are 1. next period asset level and 2. current persistent income) containing the updated bond price function
#4. qSum: a two dimensional array (axes are 1. next period asset level and 2. nexy period persistent income) containing the the value of the RHS terms
#in the functional equation describing the bond price function, i.e. qGrid[i,j]=E[1/R*qSum[i',j]|i]
#5. solveMarkH: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the highest level of the m shock has been solved during the current iteration
#6. solveMarkL: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the lowest level of the m shock has been solved during the current iteration.
#7. feasibleSolutionH: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the highest level of the m shock has been solved during the current iteration and that solution resulted in strictly positive consumption.
#8. feasibleSolutionL: a two dimensional array (axes are 1. current debt level and 2. current persistent income) containing Boolean variables indicating
#whether the government's problem at the lowest level of the m shock has been solved during the current iteration and that solution resulted in strictly positive consumption.
#9. maxAlwaysDefInd: a one dimensional array (axes are 1. current persistent income) containing the maximum asset level
#at which the government defaults under all possible realizations of the m shock

struct LongTermBondRDSaveUpdate{F<:Real}
    VF::VFuncSave{F}
    solveMarkH::Array{Bool,2}
    solveMarkL::Array{Bool,2}
    feasibleSolutionH::Array{Bool,2}
    feasibleSolutionL::Array{Bool,2}
end



struct VDFuncSave{F<:Real}
    vDInitFlow::Array{F,1}
    vDFutFlow::Array{F,1}
    vDInitGrid::Array{F,1}
    vDFutGrid::Array{F,1}
    EVDGrid::Array{F,1}
end

struct DPoliciesSave{F<:Real,S<:Integer}
    defThreshold::Array{F,2}
    firstRepayInd::Array{S,2}
    alwaysDefault::Array{Bool,2}
    neverDefault::Array{Bool,2}
    attemptedSqRootNeg::Array{Bool,2}
end


struct LongTermBondRDSaveExtension{F<:Real,S<:Integer}
    yDefGrid::Array{F,1}
    aExtPoints::S
    aGrid::Array{F,1}
    aGridIncr::Array{F,2}
    qGrid::Array{F,2}
    netRevM0A0Grid::Array{F,2}
    consM0Grid::Array{F,3}
    VFExtension::VFuncSave{F}
    VF::VFuncSave{F}
    VFD::VDFuncSave{F}
    polExtension::PoliciesSave{F,S}
    qSumExtension::Array{F,2}
    pol::PoliciesSave{F,S}
    dPol::DPoliciesSave{F,S}
    apSearchParams::GridSearchParams{S}
    noConvergeMark::Array{Bool,1}
end


struct LongTermBondRDExtensionUpdate{F<:Real,S<:Integer}
    VFExtension::VFuncSave{F}
    VF::VFuncSave{F}
    qGridExtension::Array{F,2}
    qSumExtension::Array{F,2}
    maxAlwaysDefInd::Array{S,1}
    solveMarkH::Array{Bool,2}
    solveMarkL::Array{Bool,2}
    feasibleSolutionH::Array{Bool,2}
    feasibleSolutionL::Array{Bool,2}
end




#######################################################################
# Policy function structs 
#######################################################################


struct ConsumptionPolicy{F<:Real}
    consRGrid::Array{F,2}
    consGrid::Array{F,2}
    consNoPen::Array{F,2}
    consNoPenD::Array{F,1}
    consRVar::Array{F,2}
    consVar::Array{F,2}
    consDVar::Array{F,1}
end

struct ConsumptionDivYPolicy{F<:Real}
    consRDivYGrid::Array{F,2}
    consDivYGrid::Array{F,2}
    consDivYDInitGrid::Array{F,1}
    consDivYDFutGrid::Array{F,1}
    consRDivYVar::Array{F,2}
    consDivYVar::Array{F,2}
    consDivYDInitVar::Array{F,1}
    consDivYDFutVar::Array{F,1}
end

struct UtilityPolicy{F<:Real}
    uGrid::Array{F,2}
    uRGrid::Array{F,2}
    uDInitGrid::Array{F,1}
    uDFutGrid::Array{F,1}
    uNoPen::Array{F,2}
    uNoPenD::Array{F,1}
    utYObs::Array{F,2}
    utYObsD::Array{F,1}
end

struct OutputPolicy{F<:Real}
    yObs::Array{F,2}
    yObsD::Array{F,1}
end

struct APrimeDefPolicy{F<:Real,S<:Integer}
    apRGrid::Array{F,2}
    apRDivYGrid::Array{F,2}
    defaultProb::Array{F,2}
    repayProb::Array{F,2}
    apProbability::Array{SparseMatrixCSC{F,S},1}
    lastAlwaysDefInd::Array{S,1}
end


struct MeanPolicies{F<:Real,S<:Integer}
    cPol::ConsumptionPolicy{F}
    cDivYPol::ConsumptionDivYPolicy{F}
    apDefPol::APrimeDefPolicy{F,S}
    uPol::UtilityPolicy{F}
    yPol::OutputPolicy{F}
end




#######################################################################
# Simulation structs 
#######################################################################


struct MCPaths{F<:Real,S<:Integer}
    cSim::Array{F,2}
    ySim::Array{F,2}
    mSim::Array{F,2}
    ymSim::Array{F,2}
    aSim::Array{F,2}
    apSim::Array{F,2}
    tbSim::Array{F,2}
    spreadSim::Array{F,2}
    mvSim::Array{F,2}
    defSim::Array{Bool,2}
    inDefSim::Array{Bool,2}
    yIndSim::Array{S,2}
    aIndSim::Array{S,2}
    apIndSim::Array{S,2}
    noDefDuration::Array{S,2}
    inSample::Array{Bool,2}
    inDefSample::Array{Bool,2}
end


