import numpy as np

h_fus = 2.040
b_fus = 1.19
l_fus = 3.622
kin_visc = 1.4207 * 10**-5
Cf_fus = 0.004
d_bars = 0.02
V_cruise = (100/3.6) * 1.2
rho_cruise = 1.17315

def fuselage_drag(l, b, h, rho):
    F = 2 * (l / (b + h))
    # C_D0 = (3 * F + 4.5 * (1 / F)**0.5 + 21 * (1 / F)**2) * C_f
    # S_0 = 0.25*h*b*np.pi
    # f = C_D0*S_0
    # drag_fus = f * 0.5 * rho * V_cruise**2

    C_D = (1 +1.5*(1/F)**1.5 + 7*(1/F)**3) * Cf_fus
    a1_2 = h/2
    b1_2 = b/2
    C_ellipse = np.pi * (a1_2 + b1_2) * ((3 * (a1_2-b1_2)**2)/((a1_2 + b1_2)**2 * (np.sqrt(-3*(((a1_2-b1_2)**2)/((a1_2 + b1_2)**2)) +4 ) + 10)) + 1)
    S_wet = 0.75*C_ellipse*l
    drag_fus = C_D*0.5*rho*V_cruise**2 * S_wet
    return drag_fus

def cylinder_drag_par(l, d, rho, V):  # circular cross-section is parallel to the flow
    global C_D
    F = l/d
    if F < 2:
        C_D = ((1.17-0.81)/(0-2))*F + 1.17  # DID UNIT TEST ON THIS !!!!! SLOPE WAS WRONG
    elif F > 2:
        C_D = 0.81

    S_wet = np.pi*d*l
    Re = (V*d)/kin_visc
    friction_drag = (1.328 * 0.5*V_cruise**2 * l)/(np.sqrt(Re))
    S0 = np.pi*(d/2)**2
    drag_cyl = C_D*0.5*rho*V**2 * S0 + friction_drag
    return drag_cyl

def cylinder_drag_perp(l, d, rho, V):  # circular cross-section is perpendicular to the flow
    global C_D
    Re = (V*d)/kin_visc
    if 2*10**3 < Re < 2*10**4:
        C_D = 0.4*np.log10(Re) - 0.920412
    elif 2*10**4 < Re < 2*10**5:
        C_D = 1.2
    elif 2*10**5 < Re < 5*10**5:
        C_D = -2.261647*np.log10(Re) + 13.1891
    elif 5*10**5 < Re < 10**6:
        C_D = 0.3
    elif 10**6 < Re < 3.5*10**6:
        C_D = 0.7352*np.log10(Re) - 4.111

    drag_cyl = C_D * 0.5*rho*V**2*d*l
    return drag_cyl


# Calculate drag of upper structure
prop_structure_angles = [np.radians(45), np.radians(15), np.radians(75), np.radians(45), np.radians(90)]
prop_structure_lenghts = [1.761, 1.3, 1.3, 1.3, 1.245]
drag_prop_struct = 0
drag_perp = 0
drag_par = 0
drag_hub = 0
# for i in range(len(prop_structure_lenghts)):
#     if prop_structure_angles[i] == 90:
#         drag_prop_struct += cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho_cruise, V_cruise)*4
#     else:
#         drag_par = np.cos(prop_structure_angles[i]) * cylinder_drag_par(prop_structure_lenghts[i], d_bars, rho_cruise, np.cos(prop_structure_angles[i])* V_cruise)
#         drag_perp = np.sin(prop_structure_angles[i]) * cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho_cruise, np.sin(prop_structure_angles[i]) * V_cruise)
#         drag_prop_struct += (np.sin(prop_structure_angles[i])*drag_par + np.cos(prop_structure_angles[i])*drag_perp)*4
# drag_prop_struct += 2*cylinder_drag_par(2.490, d_bars, rho_cruise, V_cruise)
# drag_prop_struct += 12*cylinder_drag_perp(0.250, 0.2, rho_cruise, V_cruise)

for i in range(len(prop_structure_lenghts)):
    if prop_structure_angles[i] == 90:
        drag_perp += cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho_cruise, V_cruise)*4
    else:
        drag_par_angle = np.cos(prop_structure_angles[i]) * cylinder_drag_par(prop_structure_lenghts[i], d_bars, rho_cruise, np.cos(prop_structure_angles[i])* V_cruise)
        drag_perp_angle = np.sin(prop_structure_angles[i]) * cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho_cruise, np.sin(prop_structure_angles[i]) * V_cruise)
        drag_par += np.sin(prop_structure_angles[i])*drag_par_angle * 4
        drag_perp += np.cos(prop_structure_angles[i])*drag_perp_angle * 4
drag_par += 2*cylinder_drag_par(2.490, d_bars, rho_cruise, V_cruise)
drag_hub += 12*cylinder_drag_perp(0.250, 0.2, rho_cruise, V_cruise)

drag_prop_struct = drag_par + drag_perp + drag_hub

# Calculate drag of fuselage
drag_fuselage = fuselage_drag(l_fus, b_fus, h_fus, rho_cruise)

# Calculate drag of landing skid
drag_skid = 0
drag_skid_vert_bar = cylinder_drag_perp(0.3, d_bars, rho_cruise, V_cruise)
drag_skid_hor_bar = cylinder_drag_par(1.5, d_bars, rho_cruise, V_cruise)
drag_skid += 2*(drag_skid_vert_bar*2 + drag_skid_hor_bar)

# Calculate total drag
total_drag = drag_prop_struct + drag_fuselage + drag_skid

# Calculate percentage of drag to weight
perc_D_W = total_drag/(1000*9.80665) * 100

print("---Drag calculations in forward flight---")
print("Upper structure drag = ", str(drag_prop_struct), "N")
print("Fuselage drag = ", str(drag_fuselage), "N")
print("Skid drag = ", str(drag_skid), "N")
print("Total drag = ", str(total_drag), "N")
print("Drag as percentage of weight = ", str(perc_D_W), "%")
print()


# print("Perpendicular to flow:", str(round(cylinder_drag_perp(1, 0.02, rho_cruise, V_cruise), 3)), "N")
# print("Parallel to flow:", str(round(cylinder_drag_par(1, 0.02, rho_cruise, V_cruise), 3)), "N")

print("---Upper structure drag breakdown---")
print("Drag parallel flow =", str(drag_par/drag_prop_struct * 100), "%")
print("Drag perpendicular flow =", str(drag_perp/drag_prop_struct * 100), "%")
print("Drag hubs =", str(drag_hub/drag_prop_struct * 100), "%")