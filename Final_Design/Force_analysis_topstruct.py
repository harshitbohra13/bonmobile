import numpy as np


class DPTrusses:
    n = 12                      # Number of rotors
    T_vec = np.zeros(n, 1)      # Initialise a vector containing the thrust at each point

    def __init__(self, xloc, yloc):
        self.xloc = xloc
        self.yloc = yloc

# Truss_A = DPTrusses(xloc, yloc)
