import numpy as np
import matplotlib.pyplot as plt

class PropSizing():
    n_rotors = 12
    R_rotor = 0.62
    rho = 1.225
    W = 4414
    A = R_rotor**2 * np.pi
    A_T = 12*A
    C_T = 0.0102
    omega = 628.318531

    def v_h(self):
        return np.sqrt(PropSizing.W/(2*PropSizing.rho*PropSizing.A_T))

    def v_h_2(self):
        return np.sqrt(PropSizing.C_T/2)*PropSizing.R_rotor*PropSizing.omega

    def T(self):
        return PropSizing.W / PropSizing.A_T

propellor = PropSizing()
print(propellor.v_h())
print(propellor.v_h_2())
print(propellor.T())

