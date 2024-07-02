import numpy as np
from scipy.optimize import minimize


'''
Class to minimize the covariance between all proteins and beta factors
proteins are M-B = measured protein - beta factor
'''
class CorrMini():

    def __init__(self, scdata):
        
        #log counts
        self._scdata = scdata
        self._scdatalog = self._scdata.applymap(np.log10)
        
        #log libsizes
        self._libSize = scdata.sum(axis=1)
        self._logLibsize = self.__calc_libsize(scdata)
        self._starting_values = self._calc_start_values()
        self._betaSum = 0.0

    #these are also loglib values, however, they sum to zero
    def _calc_start_values(self):
        #summing to zero means we need to substruct a value Y from all values x in the loglibsize vector (substration equals division in origional space)
        logSum = np.sum(self._logLibsize)
        averageValue = logSum / len(self._logLibsize)
        scaledValue = self._logLibsize - averageValue
        return(scaledValue)
        
    def __calc_libsize(self, scdata):
        
        if any(value <= 0 for value in self._libSize):
            raise ValueError("Sample with zero libsize detected! This is not possible, stop computing.")

        logLibDict = [np.log10(x) for x in self._libSize]
        return(logLibDict)
    
    # Define the system of equations
    def _equations_covariance(self, x, eta):        
        equations = []
        #for every protein set up an equation
        for protID in range(0, len(self._scdata.columns)):
            equationTmp = 0
                    #for every cell add p*b - b**2            
            std = np.std(self._scdatalog.iloc[:, protID])
                    #for cellID in range(0, self._scdata.shape[0]):
                    #    equationTmp += ((self._scdatalog.iloc[cellID, protID] * x[cellID]) - (x[cellID]**2) / std)
            equationTmp += np.cov(np.array(self._scdatalog.iloc[:, protID] - x), np.array(x))[0,1] /std    
            equations.append(eta * equationTmp)
        
        #return the whole set of equations
        return(np.array(equations))
    
    def _equations_libsize(self, x, eta):
        
        newLibSize = np.array(self._logLibsize - x)
        libSizeVar = ( (1-eta) * np.var(newLibSize) )
        
        #equations.append( -(1-eta)*np.var((x)) )
        return(libSizeVar)

    # Define the objective function
    def _objective_function(self, x, eta):
        eqs1 = self._equations_covariance(x, eta)
        eqs2 = self._equations_libsize(x, eta)

        return (np.sum(eqs1**2) + eqs2)

    def __divide_row_by_value(self, row, value):
        return (row / value)
    
    def _constraint(self, x):
        equation = 0
        for x_i in x:
            equation += x_i
        #equation -= 1
        return equation

    def solve(self, eta):

        #we assume that we divide the counts by a beta (N = M-B)
        # the first initla guess is that beta=libsize vector
        #we do not constrain the sum of betas, like this many possible solutions exist, but it does not matter as long
        # as we find any one of them...    initial_guess = np.random.rand(5)  # Assuming you have 5 parameters
        initial_beta = self._starting_values
        
        #initial_beta = np.random.rand(len(initial_beta))
        #initial_beta = np.random.uniform(low = 1, high = 1000, size = (len(initial_beta),))
        
        # Minimize the objective function
        cons = {'type':'eq', 'fun': self._constraint}
        self._result = minimize(self._objective_function, initial_beta, method='SLSQP', args=eta,  constraints=[cons])  #method='SLSQP'), BFGS
        
        #calculate the final data
        self._trueBetaFactors = [10**x for x in self._result.x]
        #self._trueBetaFactors = [10**x for x in initial_beta]
        
        #calculate the normalized matrix
        self._normalizedScData = self._scdata.apply(lambda row: self.__divide_row_by_value(row, self._trueBetaFactors[self._scdata.index.get_loc(row.name)]), axis=1)
                
    def get_logbeta_values(self):
        return(self._result.x)
    
    def get_beta_values(self):
        return(self._trueBetaFactors)
    
    def get_normalized_data(self):
        return(self._normalizedScData)