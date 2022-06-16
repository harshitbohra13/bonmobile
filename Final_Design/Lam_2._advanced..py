from cmath import sin
from pickle import TRUE
from re import T
from composites import laminated_plate
import numpy as np

#---------CHARACTERISTICS OF A SINGLE LAYER
#material characteristics
E_11 = 255# E modlus along fibres 
E_22 = 3# E modulus 
v_12 = 0.15# poisson ratio 12, actually this difference is important v_lt vs v_tl, THIS IS MAJOR POISSON'S RATIO ALONG FIBER'S DIRECTION
v_21 = 0.15# poisson ratio 21
G_12 = 2 #Shear modulus?
G_13 = 2 #Shear modulus?
G_23 = 2 #Shear modulus?

#criteria for failure in stress
max_stress_along_fibres = 1479
max_stress_tran_to_fibres = 80
max_shear_stress = 70


#----------------------------------------ABD matrix----------------------
laminaprop = (E_11, E_22, v_12 , G_12, G_13, G_23) #?
plyt = 0.15 #thickness of one layer 

#CONFIGURATION OF LAYERS
# stack = [0,90,45,-45,-45,45,90,0] #?
stack = [157.5,112.5,67.5,22.5,22.5,67.5,112.5,157.5] #?

# stack = [0, 90,0,90] #?
plate = laminated_plate(stack, plyt=plyt, laminaprop=laminaprop)
# print("ABD:",plate.ABD)

#----------------LOADS METHOD-----------------------------------------------------
LOADS = np.array([[1],[0],[0],[0],[0],[0]]) #?
# print(LOADS)
strains  = np.dot(np.linalg.inv(plate.ABD), LOADS)
# print(strains)

#-----------Calculate strains for every layer! Consider only EVEN layers!--------

#NEED TO REPAIR THE STACK LIST CHANGES INTO AN ARRAYS AND AFFECTS THE REST STRESS METHOD

# lamina_strains_in_xy = np.zeros(len(stack))
# print(len(stack))
# print(lamina_strains_in_xy)
# for i in range(0,len(stack)): 
#     #creating list of lists with x,y,xy strains for each layer
#     lamina_strains_in_xy[i] = [strains[0,0]+strains[3,0]*(plyt*(len(stack)/2-i)-plyt/2),strains[1,0]+strains[4,0]*(plyt*(len(stack)/2-i)-plyt/2), strains[2,0]+strains[5,0]*(plyt*(len(stack)/2-i)-plyt/2)]
# print(stack)

# print(lamina_strains_in_xy)

#making ana array out of the strain list
# lamina_strains_in_xy_array = np.asarray(lamina_strains_in_xy)
# print(lamina_strains_in_xy_array)
#---------------------------------------------------------------------------------

#------------------------------STRESS METHOD--------------------------
#Defining stresses
Stress_average = np.array([[600],[0],[0]])

#Defining A matrix
stiffness_matrix = plate.ABD
A_matrix = stiffness_matrix[0:3]
A_matrix = np.array([stiffness_matrix[0][0:3], stiffness_matrix[1][0:3], stiffness_matrix[2][0:3]])
# print(A_matrix)

#Definind total thickness of a composite
t_tot = plyt * len(stack)


#Find composite strains
Strains_composite = np.dot(np.linalg.inv(A_matrix), t_tot*Stress_average) 
# print(Strains_composite)


#Initialise matrix for the layer
layers_strains = np.array([[0],[0],[0]])
# print(layers_strains)

#Find/Transform strains for every layer, every column is a one layer (e_x,e_y, y_xy)
for i in range(0,len(stack)):#len(stack)

    #angle for every layer
    theta = np.deg2rad(stack[i]) #angle
 
    #matrix inputs
    T_11 = np.cos(theta)**2
    T_12 = np.sin(theta)**2
    T_13 = np.sin(theta)*np.cos(theta)
    T_21 = np.sin(theta)**2
    T_22 = np.cos(theta)**2
    T_23 = -np.sin(theta)*np.cos(theta)
    T_31 = -2*np.sin(theta)*np.cos(theta) 
    T_32 = 2*np.sin(theta)*np.cos(theta)
    T_33 = (np.cos(theta)**2 - np.sin(theta)**2)

    #rotation matrix
    strain_rotation_matrix = np.array([[T_11,T_12,T_13], [T_21,T_22,T_23], [T_31,T_32, T_33]])

    #array of strains of every layer in layer's coordinates
    new_layer_strains = np.dot(strain_rotation_matrix,Strains_composite)

    # layers_strains = np.append(layers_strains,new_layer_strains)

    layers_strains = np.hstack((layers_strains,new_layer_strains))


#delete the first initialised column of 0,0,0
#Final matrix of strains for all layers
layers_strains = layers_strains[0:3, 1:(len(stack)+1)]

#Find stress for every layer
#Inputs for strain-stress matrix
C_11 = E_11/(1-v_12*v_21)
C_12 = v_12*E_22/(1-v_12*v_21)
C_22 = E_22/(1-v_12*v_21)
C_66 = G_12

#Definign the matrix for composite layers
strain_stress_Matrix = np.array([[C_11,C_12,0],[C_12,C_22,0],[0,0,C_66]])

layers_stress = layers_strains

#calculate the stresses
for i in range (0,len(stack)):
    layers_stress[:,i] = np.dot(strain_stress_Matrix,layers_stress[:,i])
# print("Stresses", layers_stress)

#---------------FAULURE CHECK------------------
Failure = False
#Check if any layer fails
for i in range (0,len(stack)):
    layer = i
    if  layers_stress[0,i] < max_stress_along_fibres and layers_stress[1,i] < max_stress_tran_to_fibres and layers_stress[2,i] < max_shear_stress:
        Failure = False
    else:
        Failure = True 
        break
print()
print("Does the failure occur?:", Failure)
print("Maximal stress:", np.max(layers_stress), "MPa" )
print("total thickness",t_tot, "mm")
print("What layer from the top does the failure occur (0 is the 1st layer)? :", layer)
print()




