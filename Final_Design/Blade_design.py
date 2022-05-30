from turtle import width
from xml.etree.ElementTree import QName
from zipfile import ZIP_DEFLATED
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import k_means
from scipy.integrate import quad


class Bladedesign:
    
    # input parameters
    TP=0  # if 0 then for thrust if 1 then power
    B=2
    RPM=1000 
    Powerreq=10        # [kW]
    thrust_req = 3000  # [N]
    V_cruise = 100/3.6  # [m/s]
    blade_elements = 10
    disk_loading = 300     # [N/m^2]
    n_rotors = 12
    hub = 0.05    # [m]
    R = 1
    aoa = 0.5    # [rad]
    C_l = []     # list of lift coefficients for blade elements
    A = 0

    def __init__(self):
        pass

    # Basic Equations
    def rotor_area(thrust_req, disk_loading):
        A = thrust_req/disk_loading
        return A

    def propeller_radius(A):
        R=np.sqrt(A/np.pi)
        return R

    def nondimensional_radius(r, R):
        xi= r/R
        return xi

    def speed_ratio(V_freestream, omega, R):
            labda=V_freestream/(omega*R)
            return labda

    def freestream_velocity(aoa, V_cruise):
        V_freestream=np.cos(aoa)*V_cruise      # incoming flow propeller to be calculated with aoa and cruise speed
        return V_freestream

    def angular_velocity(RPM):
        omega=RPM*0.10472
        return omega

    def radial_coordinates(blade_elements, R, hub):
        dr=(R-hub)/blade_elements
        list_radialcoors=[hub+0.5*dr]
        for i in range(blade_elements-1):   
            list_radialcoors.append(list_radialcoors[i]+width)
        return list_radialcoors, dr

    # step 1 of paper
    zeta_initial=0

    # Step 2 of paper
    def prandtlequation(f):                             # Equation 18
        F=(2/np.pi)*np.arcos(np.exp(-f))
        return F
        
    def prandtlconstant(B, zheta_initial, phi_tip):     # Equation 19
        f=((B/2)*(1-zheta_initial))/np.sin(phi_tip)
        return f    
        
    def flowangle_tip(zheta_initial, labda):            # Equation 20
        phi_tip=np.arctan(labda*(1+zheta_initial/2))
        return phi_tip    
        
    def flow_angle_element(phi_tip, xi):                # Equation 21
        phi=np.arctan(np.tan(phi_tip)/xi)
        return phi    

    # Step 3 of paper

    # a list of phi is required as well as a list of F for every blade element.

    def circulation_fucntion(phi, F):
        goniometric_part= np.cos(phi) * np.sin(phi) # has to be cross product
        G=np.cross(F, goniometric_part)
        return G

    def Wcproduct(labda, G, V, R, zeta_initial, C_l, B): # Equation 16
        Wc=4*np.pi*labda*G*V*R*zeta_initial/(C_l*B)
        return Wc
    
    # Step 4 of paper
    # in this part of the design process the airfoil characteristics of the blades is analyzed.
    # For now we use the NACA4415

    epsilon = 0.05/1.35
    alpha = 12   # [deg]

    # step 5 of paper
    # iteration of this part to minimize epsilon is required

    # step 6 of paper

    def actual_interference_factor(zeta,phi, epsilon):
        a=(zeta/2)*(np.cos(phi))**2 * (1-epsilon*np.tan(phi))
        return a

    def rotational_interference_factor(x, zeta, phi, epsilon):
        a_prime=(zeta/2*x)*np.cos(phi)*np.sin(phi) *(1+ epsilon/np.tan(phi))
        return a_prime

    def relation_labda_xi(labda, xi):
        x=xi/labda
        return x

    def total_velocity(V_freestream, a, phi):
        W=V_freestream*(1+a)/np.sin(phi)
        return W
    
    # step 7 of paper
    def chord(Wc,W):
        c=Wc/W
        return c
    
    def blade_twist(phi, alpha):
        beta=phi+alpha
        return beta

    # step 8 of paper

    def numerical_integration(G,epsilon, phi, labda):
        
        def I1_prime(G1, epsilon1, phi1, xi1):
            return 4*xi1*G1*(1-epsilon1*np.tan(phi1))

        def J1_prime(G1, epsilon1, phi1, xi1):
            return 4*xi1 *G1*(1+epsilon1/np.tan(phi1))

        def I2_prime(G1, epsilon1, phi1, labda1, xi):
            return labda1*(2*G1*(1-epsilon1*np.tan(phi1)))*(1+epsilon1/np.tan(phi1))*np.sin(phi1)*np.cos(phi1)

        def J2_prime(G1, epsilon1, phi1, x1):
            return (2*x1*G1*(1-epsilon1*np.tan(phi1)))*(1-epsilon1*np.tan(phi1))*(np.cos(phi1)**2)
        
        G1 = G
        epsilon1 = epsilon
        phi1 = phi
        labda1 = labda
        I1 = quad(I1_prime, 0, 1, args=(G1, epsilon1, phi1))
        J1 = quad(J1_prime, 0, 1, args=(G1, epsilon1, phi1))
        I2 = quad(I2_prime, 0, 1, args=(G1, epsilon1, phi1, labda1))
        J2 = quad(J2_prime, 0, 1, args=(G1, epsilon1, phi1))
        return I1, I2, J1, J2

    # step 9 of paper
    # computation of displacement velocity ratio

    def displacement_velocity_ratio_throughthrust(I1, I2, T_c, J1, J2):
        zeta = (I1/2*I2)-np.sqrt((I1/2*I2)**2-T_c/I2)
        P_c = J1*zeta+J2*zeta**2
        return zeta, P_c

    def displacement_velocity_ratio_throughpower(P_c, J1, J2, I1, I2):
        zeta = -(J1/2*J2)+np.sqrt((J2/2*J2)**2+P_c/J2)
        T_c = I1*zeta-I2*zeta**2
        return zeta, T_c

    # step 10 of paper
    # if displacement velocity ration is not within 0.1%  to the initial value=> reiterate

    # step 11 of paper
    # compute propeller efficiency and other features tbd which one are of importance

    def displacement_velocity(phi, w_n):
        v_displace=w_n/np.cos(phi)
        return v_displace


def prandtlequation(blade):                             # Equation 18
    blade.F=(2/np.pi)*np.arcos(np.exp(-blade.f))
    return blade

A = 1

bonprop = Bladedesign()


def iterate(bonprop1):
    bonprop1 = prandtlequation(bonprop1)
    return(bonprop1)


for i in range(0, 100):
    bonprop = iterate(bonprop)

