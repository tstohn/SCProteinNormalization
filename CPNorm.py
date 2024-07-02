import cplex
import numpy as np
import pandas as pd

"""
CorrNorm is a class to solve a single-cell normalization problem, where one assumes that the library size scaling factor
is independant of the protein abundance: The methods solves library-size scaling factors for every individual cell, by minimizing library-size 
differences between cells (variance of lib-size), under thge constraint that the covariance of two protein-counts equals the covaraince of the normalized counts plus
variance of the library-size factor. This constraints results from the assumption that the covariance between any protein and the lib-size
factor should be close to zero!

INPUT DATA:
-------------------
@SCData:
A dataframe of single-cells times proteins. Cells are in rows (cellIDs in the rowname), and every protein is in a column.
@eta: weight for concentratin on libsize minimization / cov(prot, normFactor) minimization
    eta towards 1 results in covariances going to zero
    eta towards 0 is only library size normalisation

PUBLIC FUNCTIONS:
-------------------
@solve: solves the Normalization problem
@get_normalized_scdata: return the origional matrix with the normalized counts
@get_sclaing_factors: returns a dataframe with scaling factors for every cell, the cellIDs are in the row names
"""
class CorrOpti():
    
    #eta from 0 to 1
    def __init__(self, scdata, eta = 0.5):
        
        #raw data: index MUST contain sample names
        self._scdata = scdata
        self._scdatalog = self._scdata.applymap(np.log10)
        #list of log(libsize) in order of samples
        self._libSize = scdata.sum(axis=1)
        self._logLibsize = self.__calc_libsize(scdata)
        self._betaSum = 0.0
        #variables to infer: the norm factors per sample 
        #name: <sampleID>_beta (same order as _libsize)
        self._betaVars = []
        self._eta = eta
        
        self.__verify_format()
        
        #cplex problem
        self.__initialize_cpx()
        #create the variable names
        self.__create_variable_names()
        #set obhjective function to minimize
        self.__set_objective()
        #add constraints
        self.__set_constraints()
        
    #PRIVATE:
    def __verify_format(self):
        #check input is a dataframe
        if(not isinstance(self._scdata, pd.DataFrame)):
            print("Input file must be a dataframe. (cell*features)")
            exit(1)
        #all cells must contain numeric values
        try:
            self._scdata.applymap(pd.to_numeric)
            return True
        except ValueError:
            print("Dataframe contains non numeric values, check format (cell*features)")
            exit(1)
            
    def __calc_libsize(self, scdata):
        
        if any(value <= 0 for value in self._libSize):
            raise ValueError("Sample with zero libsize detected! This is not possible, stop computing.")

        logLibDict = [np.log10(x) for x in self._libSize]
        print(logLibDict)
        return(logLibDict)
        
    #initialize the problem
    def __initialize_cpx(self):
        self._cpx = cplex.Cplex()
        self._cpx.set_problem_type(cplex.Cplex.problem_type.QP)
        self._cpx.objective.set_sense(self._cpx.objective.sense.minimize)
        
        #self._cpx.parameters.optimalitytarget.set(2.0)

    #add all the variables
    def __create_variable_names(self):
        for sampleID in self._scdata.index:
            betaVarName = str(sampleID) + "_beta"
            self._betaVars.append(betaVarName)
            
        #every protein should be un-correlated to the beta factor
        #therefore we add an indicator for each protein weather the covariance with the beta factor is
        
        #negative or positive (to replace a quadric constraint with indicator * covariance)
        #we however need two since ind can only be [0/1]; one for fact that constraint is neg and one for pos
        #then the final linear terms are: posInd * 1.0 * cov + negInd * (-1.0) * cov
        
        self._covariances = []
        self._negCovBeyondLimitIndicator = []
        self._posCovBeyondLimitIndicator = []

        #self._negCovIndicator = []
        #self._u = []
        #self._v = []
        for protID in range(0, len(self._scdata.columns)):
            self._covariances.append('Cov_' + str(protID))
            
            self._negCovBeyondLimitIndicator.append('I_neg_' + str(protID))
            self._posCovBeyondLimitIndicator.append('I_pos_' + str(protID))

        #    self._u.append("u_" + str(protID))
        #    self._v.append("v_" + str(protID))
        self._cpx.variables.add(names= self._covariances,lb= [-cplex.infinity] * len(self._covariances), ub= [cplex.infinity] * len(self._covariances) )


            
        # NEGATIVE INDICATORS for values past the cutoff
        n_i = len(self._negCovBeyondLimitIndicator)
        self._cpx.variables.add(
            names=self._negCovBeyondLimitIndicator,
            types=[self._cpx.variables.type.binary] * n_i,
            lb=[0] * n_i,
            ub=[1] * n_i,
            obj=[self._eta] * n_i
        )
        # POSItiVe INDICATOR
        n_i = len(self._posCovBeyondLimitIndicator)
        self._cpx.variables.add(
            names=self._posCovBeyondLimitIndicator,
            types=[self._cpx.variables.type.binary] * n_i,
            lb=[0] * n_i,
            ub=[1] * n_i,
            obj=[self._eta] * n_i
        )



        #self._cpx.variables.add(names= ["negCovCutoff"],lb= [-cplex.infinity], ub= [0] )
        #self._cpx.variables.add(names= ["posCovCutoff"],lb= [0], ub= [cplex.infinity] )

        #self._cpx.variables.add(names= self._u,lb= [0] * len(self._u), ub= [cplex.infinity] * len(self._u))
        #self._cpx.variables.add(names= self._v,lb= [0] * len(self._v), ub= [cplex.infinity] * len(self._v))

        #betaLinearizer = []
        #for beta in self._betaVars:
        #    var = "beta^2_" + str(beta)
        #    betaLinearizer.append(var)
        #self._cpx.variables.add(names= betaLinearizer,lb= [0] * len(betaLinearizer), ub= [cplex.infinity] * len(betaLinearizer))
        #self._cpx.variables.add(names=["pos_threshold"], lb=[0], ub=[cplex.infinity] )
        #self._cpx.variables.add(names=["neg_threshold"], lb=[-cplex.infinity], ub=[0] )
        
        

    def __set_objective(self):
        
        
        #1. st constraint
        #___________________
        #minimizing the variance of library size means:
        #min: -2beta*log(libsize) + beta^2
        
        #add linear term: -2beta*d with d = log(libsize)
        coef = []    
          
        #the bounds are defined by differences between the two most extreme library size
        #we can MAXIMALLY decrease the highest libsize to the lowest and vice versa
        lowerBounds = [ (np.min(self._logLibsize) - np.max(self._logLibsize)) ] * len(self._betaVars)
        upperBounds = [ (np.max(self._logLibsize) - np.min(self._logLibsize)) ] * len(self._betaVars)
        for x in range(0, len(self._betaVars)):
            #self._cpx.objective.set_linear_coefficients(self._libsize[x], self._betaVars[x], -2.0)
            # we have to add the weird term of 1/2 bcs cplex for some reasons divides the quadric terms by 2
            coef.append(-2.0*self._logLibsize[x]*(1/2)*(1-self._eta))
        self._cpx.variables.add(obj = coef, names = self._betaVars, lb = lowerBounds, ub = upperBounds)

        #add quadric terms: term = +beta^2
        #minimize 
        for beta in self._betaVars:
            self._cpx.objective.set_quadratic_coefficients(beta, beta, 1.0*(1-self._eta)) 
            
        #2.) constraint
        #___________________
        #minimize also the covariance
        
        for protID in range(0, len(self._scdata.columns)):
            covtmp = "Cov_" + str(protID)
            
            coef = []
            std = np.std(self._scdatalog.iloc[:, protID])
            for cellID in range(0, self._scdata.shape[0]):
                #scale ONLY the actual protein counts
                coef.append( self._scdatalog.iloc[cellID, protID]/std * (+1.))
            #covaraince between protein & beta
            #NEGATIVE THRESHOLD
            #G means linear equ >= 0
            linearExpr = cplex.SparsePair(ind=self._betaVars, val=coef)
            self._cpx.objective.set_linear(self._betaVars, coef, (self._eta)) 
            
            quadricExpr = cplex.SparseTriple(ind1=self._betaVars, ind2=self._betaVars, val=[(-1.0)/std] * len(self._betaVars))
            
            
            
            self._cpx.quadratic_constraints.add(
                lin_expr=linearExpr,
                quad_expr=quadricExpr,
                sense="E",  
                rhs=0,
                name = "CovarianceEquation_" + str(protID)
            )
        
        
        
        
        
        
        
        # minimize the thresholds for cov(beta,protein)
        #self._cpx.objective.set_linear(  zip(   ["negCovCutoff"],  (-1) * [self._eta]    )) 
        #self._cpx.objective.set_linear(  zip(   ["posCovCutoff"],  (+1) * [self._eta]    )) 
        #self._cpx.objective.set_linear(zip(self._u, len(self._u) * [self._eta])) 
        #self._cpx.objective.set_linear(zip(self._v, len(self._v) *[self._eta]))
        #self._cpx.objective.set_linear("neg_threshold", -1.0 * (self._eta)) 


    def __set_constraints(self):
        
        #old constraint for squeezing the covariances between thresholds...
        '''
                linearExpr = None
        quadricExpr = None
        for protID in range(0, len(self._scdata.columns)):
            var = self._betaVars
            coef = []
            std = np.std(self._scdatalog.iloc[:, protID])
            for cellID in range(0, self._scdata.shape[0]):
                #scale ONLY the actual protein counts
                coef.append( self._scdatalog.iloc[cellID, protID]/std * (+1.))
            #covaraince between protein & beta
            #NEGATIVE THRESHOLD
            linearExpr = cplex.SparsePair(ind=var + ["neg_covThres"], val=coef + [-1.0])
            quadricExpr = cplex.SparseTriple(ind1=var, ind2=var, val=[(-1.0)/std] * len(var))
            self._cpx.quadratic_constraints.add(
                lin_expr=linearExpr,
                quad_expr=quadricExpr,
                sense="G",  
                rhs=0.,
                name = "cov_inverted"
            )
        '''
        
        
        #constraints:
        
        # 1.) BETA-FACTORS MUST SUM TO ZERO
        #________________________________
        
        #otherwise many solutions would be possible
        lin_exprs = [] #we have only one lin_expr still must be stored as list element here
        var = []
        coef = []
        for betaName in self._betaVars:
            var.append(betaName) #names of all the beta variables
            coef.append(1.0) #sum of beta == 0 -> coefficients are 1.0
        lin_exprs.append([var, coef])
        #set this linear equation to zero
        self._cpx.linear_constraints.add(
            lin_expr=lin_exprs,
            senses=["E"] * len(lin_exprs),
            rhs=[self._betaSum] * len(lin_exprs)
        )
        
        #2.) 
        # INDICATOR CONSTRAINT FOR NEG/POS COV(BETA, PROT)
        #________________________________

        # if cov(beta, protX) < 0 => ind == -1
        # else if cov(beta, protX) > 0 => ind == +1
        #minimizing cov(beta, prot) is:
        #min: cov(protX, beta) = sum([protX] * [beta]) + sum((-1) * beta^2)
        
        '''linearExpr = None
        quadricExpr = None

        for protID in range(0, len(self._scdata.columns)):
            covtmp = "Cov_" + str(protID)
            
            coef = []
            std = np.std(self._scdatalog.iloc[:, protID])
            for cellID in range(0, self._scdata.shape[0]):
                #scale ONLY the actual protein counts
                coef.append( self._scdatalog.iloc[cellID, protID]/std * (+1.))
            #covaraince between protein & beta
            #NEGATIVE THRESHOLD
            #G means linear equ >= 0
            linearExpr = cplex.SparsePair(ind=self._betaVars + [covtmp], val=coef + [+1.0])
            quadricExpr = cplex.SparseTriple(ind1=self._betaVars, ind2=self._betaVars, val=[(-1.0)/std] * len(self._betaVars))
            self._cpx.quadratic_constraints.add(
                lin_expr=linearExpr,
                quad_expr=quadricExpr,
                sense="E",  
                rhs=0,
                name = "CovarianceEquation_" + str(protID)
            )
            
            # POSITIVE INDICATOR: _posCovBeyondLimitIndicator // _covariances
            covTmp = self._covariances[protID]
            indTmp = self._posCovBeyondLimitIndicator[protID]
            name = str(self._posCovBeyondLimitIndicator[protID])
            constr = cplex.SparsePair(ind=[covTmp], val=[1.])
            self._cpx.indicator_constraints.add(
                indvar=indTmp,
                complemented=1,
                rhs=std, sense='L',
                lin_expr=constr,
                name=name
            )
                
            # NEGATIVE INDICATOR: _negCovBeyondLimitIndicator // _covariances
            covTmp = self._covariances[protID]
            indTmp = self._negCovBeyondLimitIndicator[protID]
            name = str(self._negCovBeyondLimitIndicator[protID])
            constr = cplex.SparsePair(ind=[covTmp], val=[1.])
            self._cpx.indicator_constraints.add(
                indvar=indTmp,
                complemented=-1,
                rhs=-std, sense='G',
                lin_expr=constr,
                name=name
            )
        '''
            
        '''
            #neg <= covarnaice
            lin_exprs = []
            var = ['negCovCutoff', covtmp]
            coef = [1.0, -1.0]
            linearExpr = cplex.SparsePair(ind=var, val=coef)
            #set this linear equation to zero
            self._cpx.linear_constraints.add(
                lin_expr=[linearExpr],
                senses=["L"],
                rhs=[0],
                names = ["NEG"]
            )
                  
            #cov <= pos      
            lin_exprs = []
            var = ['posCovCutoff', covtmp]
            coef = [1.0, -1.0]
            linearExpr = cplex.SparsePair(ind=var, val=coef)
            #set this linear equation to zero
            self._cpx.linear_constraints.add(
                lin_expr=[linearExpr],
                senses=["G"],
                rhs=[0],
                names = ["POS"]
            )
        '''
            
        '''linearExpr = None
        quadricExpr = None
        for protID in range(0, len(self._scdata.columns)):
            vProteinThres = "v_" + str(protID)
            var = self._betaVars
            coef = []
            std = np.std(self._scdatalog.iloc[:, protID])
            for cellID in range(0, self._scdata.shape[0]):
                #scale ONLY the actual protein counts
                coef.append( self._scdatalog.iloc[cellID, protID]/std * (+1.))
            #covaraince between protein & beta
            #NEGATIVE THRESHOLD
            #G means linear equ >= 0
            linearExpr = cplex.SparsePair(ind=var + [vProteinThres], val=coef + [+1.0])
            quadricExpr = cplex.SparseTriple(ind1=var, ind2=var, val=[(-1.0)/std] * len(var))
            self._cpx.quadratic_constraints.add(
                lin_expr=linearExpr,
                quad_expr=quadricExpr,
                sense="G",  
                rhs=0,
                name = "v_" + str(protID)
            )'''
        
        # POSITIVE INDICATOR
        '''n_i = len(self._posCovIndicator)
        self._cpx.variables.add(
            names=self._posCovIndicator,
            types=[self._cpx.variables.type.binary] * n_i,
            lb=[0] * n_i,
            ub=[1] * n_i
        )
        # NEGATIVE INDICATOR
        n_i = len(self._negCovIndicator)
        self._cpx.variables.add(
            names=self._negCovIndicator,
            types=[self._cpx.variables.type.binary] * n_i,
            lb=[0] * n_i,
            ub=[1] * n_i
        )
            
        # POSITIVE INDICATOR
        for num in range(0, len(self._posCovIndicator)):
            covTmp = self._covariance[num]
            indTmp = self._posCovIndicator[num]
            name = str(self._posCovIndicator[num])
            constr = cplex.SparsePair(ind=[covTmp], val=[1.])
            self._cpx.indicator_constraints.add(
                indvar=indTmp,
                complemented=1,
                rhs=0., sense='G',
                lin_expr=constr,
                name=name
            )
        # NEGATIVE INDICATOR
        for num in range(0, len(self._negCovIndicator)):
            covTmp = self._covariance[num]
            indTmp = self._negCovIndicator[num]
            name = str(self._negCovIndicator[num])
            constr = cplex.SparsePair(ind=[covTmp], val=[1.])
            self._cpx.indicator_constraints.add(
                indvar=indTmp,
                complemented=-1,
                rhs=0., sense='L',
                lin_expr=constr,
                name=name
            )'''

    def write_problem(self):
        self._cpx.write("./test.lp")
        
    def solve(self):
        self._cpx.solve()
        
        #get all the beta factors
        self._logBetaFactors = self.get_beta_values()
        #take exp of all values (we worked in log10 space)
        self._trueBetaFactors = [10**x for x in self._logBetaFactors]
        
        #calculate the normalized matrix
        self._normalizedScData = self._scdata.apply(lambda row: self.__divide_row_by_value(row, self._trueBetaFactors[self._scdata.index.get_loc(row.name)]), axis=1)

    def get_beta_values(self):
        
        assert self._cpx.solution.is_primal_feasible()
        
        variableNames = self._cpx.variables.get_names()
        betaVars = [var for var in variableNames if var.endswith('_beta')]
        
        variableValues = self._cpx.solution.get_values(betaVars)
        return(variableValues)
    
    def get_covariance_thresholds(self):
        assert self._cpx.solution.is_primal_feasible()
        
        variableNames = self._cpx.variables.get_names()
        covThresholds = [var for var in variableNames if var.endswith('_covThres')]
        
        variableValues = self._cpx.solution.get_values(covThresholds)
        return(variableValues)
        
    def __divide_row_by_value(self, row, value):
        return (row / value)

    def get_librarySize_variance(self):
        
        #divide origional library sizes by this value
        scaledLibsize = self._libSize / self._trueBetaFactors

        #then calculate var of those values
        return(np.var(scaledLibsize))
        
    #log scale
    def get_mean_of_absolute_covariances_log(self):
        absCovSum = 0
        for proteinID in range(0, self._scdata.shape[1]):
            
            std = np.std(self._scdatalog.iloc[:, proteinID])
            proteinCounts = self._scdatalog.iloc[:, proteinID]
            print("Protein counts")
            print("_____________")
            print(proteinCounts)
            covTmp = np.cov( proteinCounts, self._logBetaFactors)[0][1]/std
            print("Calced Cov: ")
            print(covTmp)
            absCovSum += np.abs(covTmp)
        return(absCovSum/self._scdata.shape[1])
    
    #my weird sclaing
    """def get_mean_of_absolute_covariances(self):
        absCovSum = 0
        for proteinID in range(0, self._scdata.shape[1]):
            proteinCounts = self._normalizedScData.iloc[:, proteinID]
            covTmp = np.cov(proteinCounts, self._trueBetaFactors)[0][1]/np.std(proteinCounts)
            absCovSum += np.abs(covTmp)
        return(absCovSum/self._scdata.shape[1])"""
    
    def get_average_corrcoef(self):
        absCovSum = 0
        for proteinID in range(0, self._scdata.shape[1]):
            proteinCounts = self._normalizedScData.iloc[:, proteinID]
            print("Protein counts")
            print("_____________")
            print(proteinCounts)
            covTmp = np.corrcoef(proteinCounts, self._trueBetaFactors)[0][1]
            print("Aclced Cov: ")
            print(covTmp)
            absCovSum += np.abs(covTmp)
        return(absCovSum/self._scdata.shape[1])
    
    def get_normalized_data(self):
        return(self._normalizedScData)
        