import numpy as np


class DPTrusses:
    n = 12                      # Number of rotors
    T_vec = np.zeros(n, 1)      # Initialise a vector containing the thrust at each point
    Lbc = 1                     # [m], length of beam BC
    Lbd = 1                     # [m], length of beam BD
    Ldc = np.sqrt(2)            # [m], length of beam DC
    theta_a = 45*np.pi/180      # [rad], angle with which beam AB is mounted to fuselage
    phi_cbd = 45*np.pi/180      # [rad], angle between diagonal CBD and BC
    phi_a = 45*np.pi/180        # [rad], angle between AB and horizontal
    T_tq = 3200                 # [N], Quarter thrust

    def __init__(self, xloc, yloc):
        self.xloc = xloc
        self.yloc = yloc

    def BCDangles(self, Lbc, Lbd, Ldc):
        theta_b = np.arccos((Lbc**2 + Lbd**2 - Ldc**2) / (2*Lbc*Lbd))
        theta_c = np.arccos((Lbc ** 2 + Lbd ** 2 - Ldc ** 2) / (2 * Lbc * Lbd))
        theta_d = np.arccos((Lbc ** 2 + Lbd ** 2 - Ldc ** 2) / (2 * Lbc * Lbd))
        return theta_b, theta_c, theta_d

    def full_quarter_forces(self):
        theta_b, theta_c, theta_d = self.BCDangles(self.Lbc, self.Lbd, self.Ldc)
        Fa_v = self.T_tq
        Fa = Fa_v/np.sin(self.phi_a)
        Fa_h = Fa*np.cos(self.phi_a)
        F_cf = -Fa_h*np.cos(self.theta_a)
        F_de = -Fa_h*np.sin(theta_b)
        return Fa, F_cf, F_de


# Truss_A = DPTrusses(xloc, yloc)
