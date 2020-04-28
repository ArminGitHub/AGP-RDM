"""
===========================================================================
A class containing all the methods to compute matrix elements of n-pair 
RDMs over AGP. The methods in this class use sumESP algorithm to compute 
the matrix elements.

For details see our paper: 
J. Chem. Phys. 151, 184103 (2019): https://doi.org/10.1063/1.5127850

--------------------------------------------------------------------------
Last modified: April 28, 2020
Author: Armin Khamoshi
===========================================================================
"""
import numpy as np

class AGP_wavefunction:

    def __init__(self, Gems, NumPairs):
        """
        The constructor of the class.

        Gems: array of geminal coefficients. Must be a numpy array of 
            size M, where M is the system's size. 
        NumPairs: Number of pairs. An integer 0 <= NumPairs <= M. 
        """
        self.eta = Gems
        self.np = NumPairs
    
    def Norm(self):
        """
        Computes the norm of AGP, <AGP|AGP>. This is computed by the 
        sumESP algorithm.
        """
        return sumESP(self.eta**2, self.np)

    def RDM(self, indices):
        return sumESP(self.eta**2, self.np, indices-1)        


"""
===========================================================================
Implementation of the sumESP algorithm. 

Let X = {x1,x2,...xM} and 0 <= N <= M. The algorithm computes the sum:

Sum[ i1 <... < iN ] X(i1)*X(i2) .... *X(iN)

with polynomial complexity.

See our paper for details:
J. Chem. Phys. 151, 184103 (2019): https://doi.org/10.1063/1.5127850

Inputs:--------------------------------------------------------------------
   X: is a numpy array of dimension M. 
   N: is the number of products in each summand. (See above.) Must be an 
    integer.
   indices (optional): Contains a list of indices of those elements in X 
    that need to be eliminated in the sum. An array of integers.

--------------------------------------------------------------------------
Last modified: April 28, 2020
Author: Armin Khamoshi
===========================================================================
"""
def sumESP(X, N, indices = None):

    # By definition
    if N == 0:
        return 1

    # If 'indices' is not empty remove the, corresponding elements 
    # from 'X'.
    if indices is not None :
        X = np.delete(X, indices)
     
    # Finds M and allocates memory for the matrix 'S':
    M = len(X)
    S = np.zeros([M, N+1])

    # The sumESF algorithm
    S[:, 0] = 1
    S[0, 1] = X[0]
    for i in range(2, M+1):
        for j in range( max(1, i+N-M), min(i, N) + 1 ):
            S[i-1, j] = S[i-2, j] + X[i-1]*S[i-2, j-1]
    
    # The final output 
    return S[M-1, N]