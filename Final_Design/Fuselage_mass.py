from configuration_and_materials import Aluminiumbadboy
from configuration_and_materials import configure
import numpy as np
import matplotlib.pyplot as plt

Total_Fuselage_mass = 0


def prelim_fuselage_thickness(M, R, sigma_ult):
    t_min = M / sigma_ult / np.pi / R**2 
    return(t_min)


#Constrains
R = 1 #m radius of fuselage
sigma_ult = Aluminiumbadboy.y_strength_mat #MPa ultimate stress
M = 5*10**6 #Nm Moment from center of fuselage
density = Aluminiumbadboy.Density #kg/m3 Density of the material

#base design
t_min = prelim_fuselage_thickness(M, R, sigma_ult) 
mass_fuselage = density * 2 * np.pi * R * t_min

#1st Iteration
def tmin(M, theta, alpha, R, sigma_ult):
    t = (M * np.sin(alpha-theta)/np.pi/R**2/sigma_ult)
    return(np.abs(t))

alpha = np.arange(0, 360, 0.1)
alpha = np.deg2rad(alpha)
theta = np.arange(-15, 15, 3)
theta = np.deg2rad(theta)
T = np.zeros((theta.size, alpha.size))

i = 0
for thet in theta:
    j = 0    
    for alp in alpha:
        T[i][j] = (tmin(M, thet, alp, R, sigma_ult))
        j += 1
    i += 1


# plt.plot(np.rad2deg(alpha), t[15])
# plt.ylabel("t [m]")
# plt.xlabel("alpha [rad]")

# plt.show()

#2nd iteration
def Ix1x1(t, a, alpha, R):
    Ixixi = t * a**3 * np.sin(alpha+np.deg2rad(90))**2/12 
    steiner = t * a * R**2 * np.sin(alpha)**2
    return(Ixixi + steiner)

def Iy1y1(t, a, alpha, R):
    Iyiyi = t * a**3 * np.cos(alpha+np.deg2rad(90))**2/12 
    steiner = t * a * R**2 * np.cos(alpha)**2
    return(Iyiyi + steiner)

def normalized_stress_t(t, theta, alpha, Ixx, Iyy, sigma_y_mat):
    sigma_y = abs(M*R*((np.cos(theta)*np.sin(alpha)/Ixx) - (np.sin(theta)*np.cos(alpha)/Iyy))/sigma_y_mat)
    tnew = t*sigma_y 
    return(sigma_y, tnew)

alpha = np.rad2deg(np.deg2rad(alpha))
beta = 50
delta_s =  R * np.deg2rad(beta)
Ixx = 0
Iyy = 0
Tnew = np.zeros((theta.size, alpha.size))
SIGMA_N = np.zeros((theta.size, alpha.size))

for i in range(alpha.size):
    Ixx += Ix1x1(T[7][i], delta_s, alpha[i], R)
    Iyy += Iy1y1(T[7][i], delta_s, alpha[i], R)

i=0
for thet in theta:
    j = 0    
    for alp in alpha:
        (SIGMA_N[i][j], Tnew[i][j]) = normalized_stress_t(T[i][j], thet, alp, Ixx, Iyy, sigma_ult)
        j += 1
    i += 1

plt.plot(np.rad2deg(alpha), Tnew[4]*1000)
plt.ylabel("t [mm]")
plt.xlabel("alpha [rad]")
plt.show()

x=1

