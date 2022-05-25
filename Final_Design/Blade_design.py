from turtle import width
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import k_means




class Bladedesign:
    
    B=2
    RPM=1000 
    Powerreq=10       #[kW]
    thrust_req=3000 # [N]
    V_cruise=100/3.6  # [m/s]
    blade_elements=10
    disk_loading=300     # [N/m^2]
    n_rotors=12
    hub= 0.05    # [m]
    R=1
    def radial_coordinates(blade_elements, R, hub):
        width=(R-hub)/blade_elements
        list_radialcoors=[hub+0.5*width]
        for i in range(blade_elements-1):   
            list_radialcoors.append(list_radialcoors[i]+width)
        return list_radialcoors
    print(radial_coordinates(blade_elements, R, hub))
    def rotor_area(thrust_req, disk_loading):
        A=thrust_req/disk_loading
        return A

    def propeller_radius(A):
        R=np.sqrt(A/np.pi)
        return R

    # def freestream_velocity():
    #     V_freestream=      #incoming flow propeller to be calculated with aoa and cruise speed
    #     return V_freestream

    def prandtlconstant(B, zheta, phi_tip):
        f=((B/2)*(1-zheta))/np.sin(phi_tip)
        return f
    
    def prandtlequation(f):
        F=(2/np.pi)*np.arcos(np.exp(-f))
        return F

    def speed_ratio(V_freestream, omega, R):
        labda=V_freestream/(omega*R)
        return labda

    def flowangle_tip(zheta, labda):
        phi_tip=np.arctan(labda*(1+zheta/2))
        return phi_tip
    def flow_angle_element(phi_tip, xi):
        phi=np.arctan(np.tan(phi_tip)/xi)
        return phi
    def nondimensional_radius(r, R):
        xi= r/R
        return xi