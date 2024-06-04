import cplex
import numpy as np
import pandas as pd

"""
CorrNorm is a class to solve a single-cell normalization problem, where one assumes that the library size scaling factor
is independant of the protein abundance: The methods solves library-size scaling factors for every individual cell, by minimizing library-size 
differences between cells (variance of lib-size), under thge constraint that the covariance of two protein-counts equals the covaraince of the normalized counts plus
variance of the library-size factor. This constraints results from the assumption that the covariance between any protein and the lib-size
factor should be close to zero!

@SCData:
A dataframe of single-cells times proteins. Cells are in rows (cellIDs in the rowname), and every protein is in a column.

@solve: solves the Normalization problem

@get_normalized_scdata: return the origional matrix with the normalized counts

@get_sclaing_factors: returns a dataframe with scaling factors for every cell, the cellIDs are in the row names
"""
class CorrNorm():
    
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

    #add all the variables
    def __create_variable_names(self):
        for sampleID in self._scdata.index:
            betaVarName = str(sampleID) + "_beta"
            self._betaVars.append(betaVarName)
        
    def __set_objective(self):
        
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
        for beta in self._betaVars:
            self._cpx.objective.set_quadratic_coefficients(beta, beta, 1.0*(1-self._eta))
            
        #additional variables to 'squeeze' the cov(prot, beta) as close asposisble to zero
        self._cpx.variables.add(names=["neg_covThres"], lb = [-cplex.infinity], ub = [0.])
        self._cpx.variables.add(names=["pos_covThres"], lb = [0.], ub = [cplex.infinity])
        
        self._cpx.objective.set_quadratic_coefficients([("neg_covThres","neg_covThres", 1.0*self._eta)])
        self._cpx.objective.set_quadratic_coefficients([("pos_covThres","pos_covThres", 1.0*self._eta)])
        
    def __set_constraints(self):
        #constraints:
        
        # 1.) beta must sum to zero
        #________________________________
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
        
        #2.) -std(prot) < cov(prot, beta) < std(prot)
        linearExpr = None
        quadricExpr = None

        for protID in range(0, len(self._scdata.columns)):
            var = self._betaVars
            coef = []
            for cellID in range(0, self._scdata.shape[0]):
                coef.append(self._scdatalog.iloc[cellID, protID] * +1.)
                
            #covaraince between protein & beta
            #NEGATIVE THRESHOLD
            linearExpr = cplex.SparsePair(ind=var + ["neg_covThres"], val=coef + [-1.0])
            quadricExpr = cplex.SparseTriple(ind1=var, ind2=var, val=[-1.0] * len(var))
            std = np.std(self._scdatalog.iloc[:, protID])
            self._cpx.quadratic_constraints.add(
                lin_expr=linearExpr,
                quad_expr=quadricExpr,
                sense="G",  
                rhs=0.,
                name = "cov_inverted"
            )
                        
            #POSITIVE THRESHOLD
            linearExpr = cplex.SparsePair(ind=var + ["pos_covThres"], val=coef + [+1.0])
            quadricExpr = cplex.SparseTriple(ind1=var, ind2=var, val=[-1.0] * len(var))
            std = np.std(self._scdatalog.iloc[:, protID])
            self._cpx.quadratic_constraints.add(
                lin_expr=linearExpr,
                quad_expr=quadricExpr,
                sense="G",  
                rhs=0.,
                name = "cov_inverted"
            )
        
    def write_problem(self):
        self._cpx.write("./test.lp")
        
    def solve(self):
        self._cpx.solve()
        
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
        

    #def get_normalized_scdata():
    
    def get_librarySize_variance(self):
        
        #get all the beta factors
        logBetaFactors = self.get_beta_values()
        
        #take exp of all values (we worked in log10 space)
        trueBetaFactors = [10**x for x in logBetaFactors]

        #divide origional library sizes by this value
        scaledLibsize = self._libSize / trueBetaFactors
        
        #then calculate var of those values
        return(scaledLibsize)
        
    #def get_mean_of_absolute_covariances():
        
        
        
    #def plot_libsizeVar_against_meanAbsCov():
        