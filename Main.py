import numpy as np
from AGP_RDM import AGP_wavefunction
%load_ext autoreload
%autoreload 2

# Define the wavefunction
eta = np.array([1,2,3,4,5])
N = 1
myAGP = AGP_wavefunction(eta, N)

# compute sumESP
myAGP.np = 2
myAGP.Norm()
myAGP.Normalize()
myAGP.Norm()

# Compute some overlap
myAGP.np = 3
L = np.array([1,4, 3])
R = np.array([4,0,0])
myAGP.RDM(L, R)
