import numpy as np
from AGP_RDM import AGP_wavefunction
%load_ext autoreload
%autoreload 2

# Define the wavefunction
eta = np.array([1,2,3,4,5])
N = 1
myAGP = AGP_wavefunction(eta, N)

# compute sumESP
myAGP.np = 4
indices = np.array([2])
myAGP.Norm()
