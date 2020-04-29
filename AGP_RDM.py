"""
===========================================================================
This files containing methods to compute the norm and matrix elements of 
n-pair RDMs over AGP. Functions therein use sumESP algorithm to 
compute the matrix elements.

For details see our paper: 
J. Chem. Phys. 151, 184103 (2019): https://doi.org/10.1063/1.5127850

--------------------------------------------------------------------------
Last modified: April 28, 2020
Author: Armin Khamoshi
===========================================================================
"""
import numpy as np

class AGP_wavefunction:

    def __init__(self, GemCoeffs, NumPairs):
        """
        The constructor of the class.

        GemCoeffs: array of geminal coefficients. Must be a numpy array 
            of size M, where M is the system's size. 
        NumPairs: Number of pairs. An integer st 0 <= NumPairs <= M. 
        """
        self.eta = GemCoeffs
        self.np = NumPairs
    
    def Norm(self):
        """
        Computes the norm of AGP, <AGP|AGP>. It is computed by the 
        sumESP algorithm.
        """
        return sumESP(self.eta**2, self.np)

    def RDM(self, L_indices, R_indices):
        """
        Computes the matrix element associated with the n-pair RDM: 
        
        <AGP|P!(p1)...P!(pn) P(q1)...P(qn)|AGP>

        Inputs:
            R_indices: The indices of pair-annihilation operators acting 
                to the right. Numpy array of length ‘n’ as defined above.
            L_indices: The indices of pair-creation operators acting 
                to the left. Numpy array of length ‘n’ as defined above.
        """

        # If 'R_indices' or 'L_indices' have repetitive elements or if their
        # lengths are not equal to each other, gives zero:
        Length = len(L_indices)
        if( Length != len( np.unique(L_indices) ) or
            Length != len( np.unique(R_indices) ) or 
            Length != len( R_indices) ):
            return 0
        else:
            new_np = self.np - Length
            Indices = np.concatenate( (L_indices, R_indices) )

        # Computes the prefactor:
        prefactor = np.prod( self.eta[Indices] )

        return prefactor*sumESP(self.eta**2, new_np, Indices)

    def Normalize(self):   
        """
        Normalizes the geminal coefficients. This is done by:
        
        eta(p) --> eta(p) / (<AGP|AGP>)^(1/(2N))
        """

        # Computes the <AGP|AGP> and scales the geminal coefficients
        AGP_norm = self.Norm()
        self.eta = self.eta * (AGP_norm)**( -1/(2*self.np) )

        # Print a message once normalization is completed.
        print('The geminal coefficients are now normalized.' )


"""
===========================================================================
Implementation of the sumESP algorithm. 

Let X = {x1,x2,...xM} and 0 <= N <= M. The algorithm computes the following 
sum:

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

    # By definition...
    if N == 0:
        return 1
    elif(N < 0 or N > len(X)):
        return 0

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