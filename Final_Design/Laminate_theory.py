from cmath import sin
from re import T
from composites import laminated_plate
import numpy as np

#material characteristics
E_11 = 255# E modlus along fibres
E_22 = 3# E modulus 
v_12 = 0.15# poisson ratio 12, actually this difference is important v_lt vs v_tl, THIS IS MAJOR POISSON'S RATIO ALONG FIBER'S DIRECTION
v_21 = 0.15# poisson ratio 21
G_12 = 2 #Shear modulus?
G_13 = 2 #Shear modulus?
G_23 = 2 #Shear modulus?

#Inputs for strain-stress matrix
C_11 = E_11/(1-v_12*v_21)
C_12 = v_12*E_22/(1-v_12*v_21)
C_22 = E_22/(1-v_12*v_21)
C_66 = G_12

strain_stress_Matrix = np.array([[C_11,C_12,0],[C_12,C_22,0],[0,0,C_66]])
#Inputs for strain transformation
theta = np.deg2rad(90) #angle
T_11 = np.cos(theta)**2
T_12 = np.sin(theta)**2
T_13 = np.sin(theta)*np.cos(theta)
T_21 = np.sin(theta)**2
T_22 = np.cos(theta)**2
T_23 = -np.sin(theta)*np.cos(theta)
T_31 = -2*np.sin(theta)*np.cos(theta) 
T_32 = 2*np.sin(theta)*np.cos(theta)
T_33 = (np.cos(theta)**2 - np.sin(theta)**2)

strain_rotation_matrix = np.array([[T_11,T_12,T_13], [T_21,T_22,T_23], [T_31,T_32, T_33]])


laminaprop = (E_11, E_22, v_12 , G_12, G_13, G_23) #?
plyt = 0.1 #thickness of one layer ?
stack = [0, 90, 90, 0 ] #?
plate = laminated_plate(stack, plyt=plyt, laminaprop=laminaprop)
# print(plate.ABD)


LOADS = np.array([[1],[1],[1],[0],[0],[0]]) #?
# print(LOADS)
strains  = np.dot(np.linalg.inv(plate.ABD), LOADS)
print(strains)

#-----------Calculate strains for every layer! Consider only EVEN layers!--------
lamina_strains_in_xy = stack
for i in range(0,len(stack)): 
    #creating list of lists with x,y,xy strains for each layer
    lamina_strains_in_xy[i] = [strains[0,0]+strains[3,0]*(plyt*(len(stack)/2-i)-plyt/2),strains[1,0]+strains[4,0]*(plyt*(len(stack)/2-i)-plyt/2), strains[2,0]+strains[5,0]*(plyt*(len(stack)/2-i)-plyt/2)]

print(lamina_strains_in_xy)

#making ana array out of the strain list
lamina_strains_in_xy_array = np.asarray(lamina_strains_in_xy)
print(lamina_strains_in_xy_array)


