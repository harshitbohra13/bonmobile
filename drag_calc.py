from cmath import pi
from re import X
import numpy as np


#inPlaneH from rotor induced
#Parasite drag
#T*alpha = H + Dp
#Symmetric, circular fuselage


#Parameters
rho = 1.225
V   = 27.77
solidity = 0.08

#Fuesalage literature
Cd0_fuesalage    = 0.5
Cd0_structures   = 0

#assuming angle of attack at -5deg
aoa = np.deg2rad(5)

#flat plate
Cd0_blade = 0.001
 
 
 

#functions for rotor drag
#C_d0 = T(permotor)/(0.5*rho*V^2*pi*r^2)
def rotor_drag(T, r, n):
    Cl = T/(n*0.5*rho*(V**2)*pi*(r**2))
    Cd0_rotor = Cd0_blade*3 + Cl*aoa   
    return(Cd0_rotor)


def rotor_drag90(T, r):
    Cdrotor = 0.05 #flat structural drag
    Cd0_rotor = Cd0_blade*3 + Cdrotor
    return(Cd0_rotor)


#functions for cd0
def Cd0_design1(rotors_number, propeller_radius, total_T):
    number_rotors = rotors_number
    radius_rotor = propeller_radius
    T = total_T
    Cd0_structures = 0.1
    Cd0_rotors = 0
    
    for i in range(0, number_rotors):
        Cd0_rotors = Cd0_rotors + rotor_drag(T/number_rotors, radius_rotor, number_rotors)   
    return(Cd0_rotors+Cd0_structures+Cd0_fuesalage)


def Cd0_design2(rotors_number, propeller_radius, total_T, hoverthrust, horprop_d):
    number_rotors = rotors_number - 1
    number_rotors90 = 1
    radius_rotor = propeller_radius
    T = total_T
    Cd0_structures = 0.1
    Cd0_rotors = 0 
    
    # for i in range(0, number_rotors):
    #     Cd0_rotors = Cd0_rotors + rotor_drag(T/number_rotors, radius_rotor, number_rotors) 
    Cd0_rotors += rotor_drag90(hoverthrust, horprop_d/2)
    return(Cd0_rotors+Cd0_structures+Cd0_fuesalage)

    
def Cd0_design3(rotors_number, propeller_radius, total_T):
    number_rotors = rotors_number
    radius_rotor = propeller_radius
    T = total_T
    Cd0_structures = 0.13
    Cd0_rotors = 0
    for i in range(0, number_rotors):
        Cd0_rotors = Cd0_rotors + rotor_drag(T/number_rotors, radius_rotor, number_rotors) 
    return(Cd0_rotors+Cd0_structures+Cd0_fuesalage)


def Cd0_design4(rotors_number, propeller_radius, total_T, hoverthrust, horprop_d):
    number_rotors = rotors_number-2
    number_rotors90 = 2
    radius_rotor = propeller_radius
    T = total_T
    Cd0_structures = 0.1
    Cd0_rotors = 0
    
    for i in range(0, number_rotors):
        Cd0_rotors = Cd0_rotors + rotor_drag(T/number_rotors, radius_rotor, number_rotors) 
    for i in range(0, number_rotors90):
        Cd0_rotors += rotor_drag90(hoverthrust, horprop_d/2)
        
    return(Cd0_rotors + Cd0_structures + Cd0_fuesalage)

def Cd0_design5(rotors_number, propeller_radius, total_T):
    number_rotors = rotors_number
    number_rotors90 = 0
    radius_rotor = propeller_radius
    T = total_T
    Cd0_structures = 0.2
    Cd0_rotors =  rotor_drag(T/number_rotors, radius_rotor, number_rotors)
    return(Cd0_rotors+Cd0_structures+Cd0_fuesalage)
    

