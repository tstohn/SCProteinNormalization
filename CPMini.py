import numpy as np
from scipy.optimize import minimize
import scipy.optimize
import pandas as pd
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy.stats import spearmanr

'''
Class to minimize the covariance between all proteins and beta factors
proteins are M-B = measured protein - beta factor
'''
class CorrMini():

    def __init__(self, scdata, startBetaFactors = None):
        
        #log counts
        self._scdata = scdata
        self._scdatalog = self._scdata.applymap(np.log10)
        
        #log libsizes
        self._libSize = scdata.sum(axis=1)

        self._logLibsize = self.__calc_libsize(scdata)
        if(startBetaFactors is None):
            self._starting_values = self._calc_start_values()
        else:
            self._starting_values = self._set_starting_values(startBetaFactors)
        self._betaSum = 0.0
        self._covSum = 0.0
        
        self._scdatalogLibnormed = self._normalize_data()
        print(self._scdatalogLibnormed)
        print(self._scdatalog)
        print(self._logLibsize)

    def _normalize_data(self):
        normedData = self._scdatalog.sub(self._logLibsize, axis=0)
        return(normedData)
        
    #we use an input dataFrame (sample_id, betaFactors) to set the starting values
    #with this we can, e.g., use TMM normalized factors to start with
    #the faction has to map sample_id rowname sof the data table to the beta factores and order them in the right way
    def _set_starting_values(self, startBetaFactors):
                
        sample_id_order = self._scdata.index.values.tolist()
        startBetaFactors['sample_id'] = pd.Categorical(startBetaFactors['sample_id'], categories=sample_id_order, ordered=True)
        startBetaFactorsSorted = startBetaFactors.sort_values('sample_id')

        startingFactors = [np.log10(x) for x in startBetaFactorsSorted['betaFactor']]
        
        #sclae sum to zero
        logSum = np.sum(startingFactors)
        averageValue = logSum / len(startingFactors)
        scaledValue = startingFactors - averageValue
                
        return(scaledValue)
    
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
    def _equations_covariance(self, x):        
        equations = []
        #for every protein set up an equation
        for protID in range(0, len(self._scdata.columns)):
            for protJD in range(0, len(self._scdata.columns)):
                if(protJD == protID): continue

                equationValue = np.corrcoef(np.array(self._scdatalog.iloc[:, protID] - x), x)[0,1]
                if(np.isnan(equationValue)):
                    print("WARNING: covariance = NAN for" + str(protID) + " " + str(np.corrcoef(np.array(self._scdatalog.iloc[:, protID] - x), np.array(x))[0,1]) + "\n")
                    equationValue = 0
                
                equations.append(equationValue)

        #return the whole set of equations
        return(np.array(equations))
    
    def _equations_libsize(self, x):
        
        newLibSize = np.array(self._logLibsize - x)
        libSizeVar = np.var(newLibSize)

        #check for nan
        if(np.isnan(libSizeVar)):
            print("WARNING: library size variance = NAN\n")
            return(0)

        #equations.append( -(1-eta)*np.var((x)) )
        return(libSizeVar)

    # Define the objective function
    def _objective_function(self, x, eta=0.8):
        
        eqs1 = self._equations_covariance(x)
        eqs2 = self._equations_libsize(x)
        
        return (eta * np.sum(eqs1**2) + (1-eta) * eqs2)

    def __divide_row_by_value(self, row, value):
        return (row / value)
    
    def _constraint_beta_sum(self, x):
        equation = 0
        for x_i in x:
            equation += x_i
        #equation -= 1
        return equation
    
    def _calculate_cell_similarity(self):
        #get the 10 closest cells for every cell after libsize normalization
        #we calculate similarity by comparing THE LOG COUNTS
        #this way we account for protein abundance differences, while still slightly weighting more abudnant proteins slightly more
        self._adjMatrix = {}
        
        # Calculate distances between each pair of rows
        for i in range(len(self._scdatalogLibnormed)):
            # Calculate Manhattan distances from row i to all other rows
            distances = self._scdatalogLibnormed.apply(lambda row: np.sum(np.abs(row - self._scdatalogLibnormed.iloc[i])), axis=1)
            # Sort indices based on distances, excluding the distance to itself
            sorted_indices = distances.argsort()
            sorted_indices = sorted_indices[sorted_indices != i]
                        
            # Store the 10 closest row indices
            self._adjMatrix[i] = sorted_indices[:self._neighbors].tolist()
    
    def _constraint_beta_similarity(self, x):
        constraints = []

        #enforce similarity of beta factors for all local neighborhoods
        for cell_i in range(len(self._starting_values)):
            adjacentCells = self._adjMatrix[cell_i]
            
            adjacentBetas = x[adjacentCells]
            adjacentLibsizeBetas = self._starting_values[adjacentCells]

            betaVariance = np.var(adjacentBetas - adjacentLibsizeBetas)
            constraints.append(betaVariance)

        return(np.array(constraints))

    def print_terms(self):
        
        print("######## START #########")
        print("Covariance Terms: ")
        eqs = self._equations_covariance(np.random.rand(len(self._starting_values)), 1.0)
        #print(eqs)
        print(str(np.sum(eqs**2)))
        
        print("Libsize Terms: ")
        eqs = self._equations_libsize(np.random.rand(len(self._starting_values)), 0.0)
        print(str(eqs))
        
        print("Difference:")
        print(np.array(self._logLibsize - self._starting_values))
        
        print("Initial Values:")
        print(self._starting_values)
        print("######## STOP #########\n\n")

    def solve(self, eta, neighbors = 50, minVariation = 0.0001):

        self._neighbors = neighbors
        self._calculate_cell_similarity()
                
        #we assume that we divide the counts by a beta (N = M-B)
        # the first initla guess is that beta=libsize vector
        #we do not constrain the sum of betas, like this many possible solutions exist, but it does not matter as long
        # as we find any one of them...    initial_guess = np.random.rand(5)  # Assuming you have 5 parameters
        initial_beta = self._starting_values
                
        #initial_beta = np.random.rand(len(initial_beta))
        #initial_beta = np.zeros(len(initial_beta))
        #initial_beta = np.random.choice([-1, 1], size=len(initial_beta))
        #initial_beta = np.random.uniform(low = 1, high = 1000, size = (len(initial_beta),))
        
        #min = np.min(self._logLibsize)
        #max = np.max(self._logLibsize)
        #diff = (max - min)/2
        #initial_beta_scatter = np.random.choice([-diff, diff], size=len(initial_beta))
        #initial_beta = initial_beta + initial_beta_scatter
        
        # Minimize the objective function
       # cons = [{'type':'eq', 'fun': self._constraint_beta_sum},
       #         {'type':'eq', 'fun': self._constraint_beta_similarity}]
        
        A = np.ones((1, len(initial_beta)))  # Row vector of ones
        b = 0  # We want the sum to be zero
        # LinearConstraint requires lb and ub to be set
        cons1 = LinearConstraint(A, lb=b, ub=b)
        cons2 = NonlinearConstraint(self._constraint_beta_similarity, 0, minVariation)
        
        if(self._neighbors > 0):
            cons = [cons1, cons2]
        else:
            cons = [cons1]

        self._result = minimize(self._objective_function, initial_beta, method='SLSQP', args=eta,  constraints=cons)  #method='SLSQP'), BFGS
        
        print("ASSERT CONSTRAINTS HOLD TRUE:")
        print(sum(self._result.x))
        #min = np.min(self._logLibsize)
        #max = np.max(self._logLibsize)
        #diff = (max - min)
        #bounds = [(-diff, diff), ]*len(initial_beta)
        
        #bounds = [(val - 1 * np.abs(diff), val + 1 * np.abs(diff)) for val in initial_beta]
        #self._result = scipy.optimize.shgo(self._objective_function, bounds,  constraints=[cons], sampling_method = 'sobol', n =100)
        #print(self._result)        

        #calculate the final data
        self._trueBetaFactors = [10**x for x in self._result.x]
        #self._trueBetaFactors = [10**x for x in initial_beta]
        
        trueStartingPoint = [10**x for x in self._starting_values]
        self._trueBetaFactorLibsizeDiff = [a - b for a, b in zip(self._trueBetaFactors, trueStartingPoint)]
        
        #calculate the normalized matrix
        self._normalizedScData = self._scdata.apply(lambda row: self.__divide_row_by_value(row, self._trueBetaFactors[self._scdata.index.get_loc(row.name)]), axis=1)
                
                
        #CHECK CELL SIMILARTIY
        for cell_i in range(len(self._starting_values)):
            #print("_______")

            #print(self._scdata.index[cell_i])
            adjacentCells = self._adjMatrix[cell_i]
            #print(adjacentCells)
            adjacentBetas = self._result.x[adjacentCells]
            libBetas = self._starting_values[adjacentCells]
            #print((adjacentBetas[0]-libBetas[0]))
                
    def get_logbeta_values(self):
        return(self._result.x)
    
    def get_beta_values(self):
        betaFrame = pd.DataFrame({"sample_id": self._scdata.index, "betaFactor" : self._trueBetaFactors})
        return(betaFrame)
    
    def get_betaLibDiff_values(self):
        betaFrame = pd.DataFrame({"sample_id": self._scdata.index, "betaFactorDiff" : self._trueBetaFactorLibsizeDiff})
        return(betaFrame)
    
    def get_normalized_data(self):
        return(self._normalizedScData)
    
    def get_librarySize_variance(self):
        
        #divide origional library sizes by this value
        print(self._libSize)
        print(self._trueBetaFactors)

        scaledLibsize = self._libSize / self._trueBetaFactors
        print(scaledLibsize)

        #then calculate var of those values
        return(np.var(scaledLibsize))
    
    def get_mean_of_absolute_covariances_log(self):
        avgCov = 0
        for proteinID in range(0, self._scdata.shape[1]):
            proteinCounts = self._normalizedScData.iloc[:, proteinID]
            std = np.std(self._normalizedScData.iloc[:, proteinID])

            covTmp = np.cov( proteinCounts, self._trueBetaFactors)[0, 1] / std
            avgCov += np.abs(covTmp)
            
        return(avgCov/self._scdata.shape[1])