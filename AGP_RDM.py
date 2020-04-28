"""
A class containing all the methods to compute matrix elements of n-pair 
RDMs over AGP. The methods in this class use sumESP algorithm to compute 
the matrix elements.

For details see our paper: 
J. Chem. Phys. 151, 184103 (2019): https://doi.org/10.1063/1.5127850

------------------------------------------------------------------------   
Last modified: April 28, 2020
Author: Armin Khamoshi
------------------------------------------------------------------------   
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
        return sumESP(self.eta**2, self.np)

    def RDM(self, indices):
        return sumESP(self.eta**2, self.np, indices)        

def sumESP(X, N, indices = None):

    # By definition
    if N == 0:
        return 1

    # If 'indices' is not empty remove the, corresponding elements 
    # from 'X'.
    if indices:
        X = np.delete(X, indices)
     
    # Finds M and allocates memory for the matrix 'S':
    M = len(X)
    S = np.zeros([M, N+1])

    # The summation algorithm 
    S[:, 0] = 1
    S[0, 1] = X[0]
    for i in range(2, M +1):
        for j in range( max(1, i+N-M), min(i, N) +1 ):
            S[i -1, j] = S[i-1 -1, j] + X[i -1]*S[i-1 -1, j-1]
    
    return S[M-1, N]