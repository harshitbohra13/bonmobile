import numpy as np
import matplotlib.pyplot as plt

m_vehicle = 1400
g = 9.80665
h_fus = 2.550
b_fus = 1.9
l_fus = 4.124
din_visc = 1.825 * 10**-5
Cf_fus = 0.004
d_bars = 0.02
V_cruise = (100/3.6) * 1.2
rho_cruise = 1.17315
l_hub = 0.250
d_hub = 0.220

def skin_coefficient(Re):
    Cf = 0.074/Re**(0.2)
    return Cf

def Reynolds(L, V):
    din_visc = 1.825 * 10 ** -5
    Re = V*L*rho_cruise/din_visc
    return Re

def fuselage_drag(l, b, h, rho):
    Re = Reynolds(l, V_cruise)
    Cf = skin_coefficient(Re)
    F = 2 * (l / (b + h))
    C_D = (1 +1.5*(1/F)**1.5 + 7*(1/F)**3) * Cf
    a1_2 = h/2
    b1_2 = b/2
    C_ellipse = np.pi * (a1_2 + b1_2) * ((3 * (a1_2-b1_2)**2)/((a1_2 + b1_2)**2 * (np.sqrt(-3*(((a1_2-b1_2)**2)/((a1_2 + b1_2)**2)) +4 ) + 10)) + 1)
    S_wet = 0.75*C_ellipse*l
    drag_fus = C_D*0.5*rho*V_cruise**2 * S_wet
    return drag_fus

def fuselage_drag_frontalArea(l, b, h, rho):
    Re = Reynolds(l, V_cruise)
    Cf = skin_coefficient(Re)
    F = 2 * (l / (b + h))
    C_D0 = (3 * F + 4.5 * (1 / F)**0.5 + 21 * (1 / F)**2) * Cf
    S_0 = 0.25*h*b*np.pi
    f = C_D0*S_0
    drag_fus = f * 0.5 * rho * V_cruise**2
    return drag_fus

def cylinder_drag_par(l, d, rho, V):  # circular cross-section is parallel to the flow
    global C_D
    F = l/d
    if F < 2:
        C_D = ((1.17-0.81)/(0-2))*F + 1.17  # DID UNIT TEST ON THIS !!!!! SLOPE WAS WRONG
    elif F > 2:
        C_D = 0.81

    Re = Reynolds(d, V)
    friction_drag = (1.328 * 0.5*V_cruise**2 * l * np.pi * d)/(np.sqrt(Re))
    S0 = np.pi*(d/2)**2
    drag_cyl = C_D*0.5*rho*V**2 * S0 + friction_drag
    return drag_cyl

def Cd_Re(Re):
    global C_D
    if Re < 2*10**3:
        C_D = 10**(-0.284*np.log10(Re) + np.log10(7.7984))
    elif 2*10**3 <= Re < 2*10**4:
        C_D = 0.3*np.log10(Re) - 0.09031
    elif 2*10**4 <= Re < 2*10**5:
        C_D = 1.2
    elif 2*10**5 <= Re < 5*10**5:
        C_D = -2.261647*np.log10(Re) + 13.1891
    elif 5*10**5 <= Re < 10**6:
        C_D = 0.3
    elif 10**6 <= Re <= 3.5*10**6:
        C_D = 0.7352*np.log10(Re) - 4.111
    elif Re > 3.5*10**6:
        C_D = 0.7
    return C_D

def cylinder_drag_perp(l, d, rho, V):  # circular cross-section is perpendicular to the flow
    global C_D
    Re = Reynolds(d, V)
    C_D = Cd_Re(Re)

    drag_cyl = C_D * 0.5*rho*V**2*d*l
    return drag_cyl

def total_drag_calc(rho, V, induced_factor):
    # Calculate drag of upper structure
    prop_structure_angles = [np.radians(45), np.radians(15), np.radians(75), np.radians(45), np.radians(90)]
    prop_structure_lenghts = [1.761, 1.3, 1.3, 1.3, 1.245]
    drag_prop_struct = 0
    drag_perp = 0
    drag_par = 0
    drag_hub = 0

    for i in range(len(prop_structure_lenghts)):
        if prop_structure_angles[i] == 90:
            drag_perp += cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho, V*induced_factor) * 4
        else:
            drag_par_angle = np.cos(prop_structure_angles[i]) * cylinder_drag_par(prop_structure_lenghts[i], d_bars,
                                                                                  rho, np.cos(
                    prop_structure_angles[i]) * V*induced_factor)
            drag_perp_angle = np.sin(prop_structure_angles[i]) * cylinder_drag_perp(prop_structure_lenghts[i], d_bars,
                                                                                    rho, np.sin(
                    prop_structure_angles[i]) * V*induced_factor)
            drag_par += np.sin(prop_structure_angles[i]) * drag_par_angle * 4
            drag_perp += np.cos(prop_structure_angles[i]) * drag_perp_angle * 4
    drag_par += 2 * cylinder_drag_par(2.490, d_bars, rho, V*induced_factor)
    drag_hub += 12 * cylinder_drag_perp(l_hub, d_hub, rho, V*induced_factor)

    drag_prop_struct = drag_par + drag_perp + drag_hub

    # Calculate drag of fuselage
    drag_fuselage = fuselage_drag(l_fus, b_fus, h_fus, rho)

    # Calculate drag of landing skid
    drag_skid = 0
    drag_skid_vert_bar = cylinder_drag_perp(0.5, d_bars, rho, V*induced_factor)
    drag_skid_hor_bar = cylinder_drag_par(1.4, d_bars, rho, V*induced_factor)
    drag_skid += 2 * (drag_skid_vert_bar * 2 + drag_skid_hor_bar)

    # Calculate total drag
    total_drag = drag_prop_struct + drag_fuselage + drag_skid
    return total_drag

def statistical_heli_drag(mass, k, rho, V):
    D_actual = total_drag_calc(rho_cruise, V, 1) / (0.5*rho*V**2)
    weight_lbs = mass *  2.205

    weight_list = np.arange(0, 8000, 100)
    k_list = [9, 2.5, 1.6, 1.4]
    k_correspondence = ["Old helicopters", "Low-drag helicopters", "Tilt-rotors", "Turboprop"]
    f_list = [[], [], [], []]
    for i in range(len(k_list)):
        for j in range(len(weight_list)):
            f = 0.092903 * k_list[i] * (weight_list[j] / 1000) ** (2 / 3)
            f_list[i].append(f)

    plt.figure(1)
    for i in range(len(k_list)):
        plt.plot(weight_list/2.205, f_list[i], label=f"k = {k_list[i]} ({k_correspondence[i]})")
    plt.plot(weight_lbs/2.205, D_actual, 'x', label="Actual mass/drag value")
    plt.xlabel("$m$ [kg]")
    plt.ylabel(f"$f$ = $D$/$q$ [$m^2$]")
    plt.legend()
    plt.show()

statistical_heli_drag(1400, 2.5, rho_cruise, V_cruise)

# Calculate drag of upper structure
prop_structure_angles = [np.radians(45), np.radians(15), np.radians(75), np.radians(45), np.radians(90)]
prop_structure_lenghts = [1.761, 1.3, 1.3, 1.3, 1.245]
drag_prop_struct = 0
drag_perp = 0
drag_par = 0
drag_hub = 0

for i in range(len(prop_structure_lenghts)):
    if prop_structure_angles[i] == 90:
        drag_perp += cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho_cruise, V_cruise)*4
    else:
        drag_par_angle = np.cos(prop_structure_angles[i]) * cylinder_drag_par(prop_structure_lenghts[i], d_bars, rho_cruise, np.cos(prop_structure_angles[i])* V_cruise)
        drag_perp_angle = np.sin(prop_structure_angles[i]) * cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho_cruise, np.sin(prop_structure_angles[i]) * V_cruise)
        drag_par += np.sin(prop_structure_angles[i])*drag_par_angle * 4
        drag_perp += np.cos(prop_structure_angles[i])*drag_perp_angle * 4
drag_par += 2*cylinder_drag_par(2.490, d_bars, rho_cruise, V_cruise)
drag_hub += 12*cylinder_drag_perp(l_hub, d_hub, rho_cruise, V_cruise)

drag_prop_struct = drag_par + drag_perp + drag_hub

# Calculate drag of fuselage
drag_fuselage = fuselage_drag(l_fus, b_fus, h_fus, rho_cruise)

# Calculate drag of landing skid
drag_skid = 0
drag_skid_vert_bar = cylinder_drag_perp(0.5, d_bars, rho_cruise, V_cruise)
drag_skid_hor_bar = cylinder_drag_par(1.4, d_bars, rho_cruise, V_cruise)
drag_skid += 2*(drag_skid_vert_bar*2 + drag_skid_hor_bar)

# Calculate total drag
total_drag = drag_prop_struct + drag_fuselage + drag_skid

# Calculate percentage of drag to weight
perc_D_W = total_drag/(m_vehicle*9.80665) * 100

print("---Drag calculations in forward flight---")
print("Upper structure drag = ", str(round(drag_prop_struct, 3)), "N    Percentage of total =", str(round(drag_prop_struct/total_drag * 100, 2)), "%")
print("Fuselage drag = ", str(round(drag_fuselage, 3)), "N            Percentage of total =", str(round(drag_fuselage/total_drag * 100, 2)), "%")
print("Skid drag = ", str(round(drag_skid, 3)), "N                 Percentage of total =", str(round(drag_skid/total_drag * 100, 2)), "%")
print("Total drag = ", str(round(total_drag, 3)), "N")
print("Drag as percentage of weight = ", str(perc_D_W), "%")
print()


fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
drag_components = ["Upper structure", "Fuselage", "Skid"]
drag_components_values = [round(drag_prop_struct, 3), round(drag_fuselage, 3), round(drag_skid, 3)]

def func(pct, allvals):
    absolute = int(np.round(pct/100.*np.sum(allvals)))
    return "{:.1f}%\n({:d} N)".format(pct, absolute)

wedges, texts, autotexts = ax.pie(drag_components_values, shadow=True, startangle=90, autopct=lambda pct: func(pct, drag_components_values),
                                  textprops=dict(color="w"))
ax.legend(wedges, drag_components,
          title="Drag Components",
          loc="center left",
          bbox_to_anchor=(0.9, 0, 0.5, 1))
plt.setp(autotexts, size=8, weight="bold")
ax.set_title(f"Total Drag = {round(total_drag, 2)} N")



# print("Perpendicular to flow:", str(round(cylinder_drag_perp(1, 0.02, rho_cruise, V_cruise), 3)), "N")
# print("Parallel to flow:", str(round(cylinder_drag_par(1, 0.02, rho_cruise, V_cruise), 3)), "N")

print("---Upper structure drag breakdown---")
print("Drag parallel flow =", str(round(drag_par,3)), "N         Percentage of total =", str(round(drag_par/drag_prop_struct * 100, 2)), "%")
print("Drag perpendicular flow =", str(round(drag_perp,3)), "N    Percentage of total =", str(round(drag_perp/drag_prop_struct * 100, 2)), "%")
print("Drag hubs =", str(round(drag_hub,3)), "N                 Percentage of total =", str(round(drag_hub/drag_prop_struct * 100, 2)), "%")
print()

fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
drag_components = ["Parallel struts", "Perpendicular struts", "Hubs"]
drag_components_values = [round(drag_par,3), round(drag_perp,3), round(drag_hub,3)]

colors = ['#7DAEFF','#003DEC','#00B0D3']

wedges, texts, autotexts = ax.pie(drag_components_values, shadow=True, startangle=90, colors=colors, autopct=lambda pct: func(pct, drag_components_values),
                                  textprops=dict(color="w"))
ax.legend(wedges, drag_components,
          title="Drag Components",
          loc="center left",
          bbox_to_anchor=(0.9, 0, 0.5, 1))
plt.setp(autotexts, size=8, weight="bold")
ax.set_title(f"Total Drag Upper Structure = {round(drag_prop_struct, 2)} N")

ld_list = np.arange(0.5, 25, 0.1)
drag_wetted_list = []
drag_frontal_list = []
drag_perp_list = []
l_hub_ex = 0.25

a1_2 = h_fus / 2
b1_2 = b_fus / 2
for i in range(len(ld_list)):
    l_value = 0.5 * ld_list[i] * (b_fus + h_fus)
    Re_value = Reynolds(l_value, V_cruise)
    Cf_value = skin_coefficient(Re_value)
    drag_wetted_list.append((1 +1.5*(1/ld_list[i])**1.5 + 7*(1/ld_list[i])**3) * Cf_value * 0.75* l_value * np.pi * (a1_2 + b1_2) * ((3 * (a1_2-b1_2)**2)/((a1_2 + b1_2)**2 * (np.sqrt(-3*(((a1_2-b1_2)**2)/((a1_2 + b1_2)**2)) +4 ) + 10)) + 1) * 0.5 * rho_cruise * V_cruise**2)
    drag_frontal_list.append((3 * ld_list[i] + 4.5 * (1 / ld_list[i])**0.5 + 21 * (1 / ld_list[i])**2) * Cf_value * 0.25*h_fus*b_fus*np.pi * 0.5 * rho_cruise * V_cruise**2)
    d_hub_ex = l_hub_ex/ld_list[i]
    drag_perp_list.append(cylinder_drag_perp(l_hub_ex, d_hub_ex, rho_cruise, V_cruise))

plt.figure(3)
plt.plot(ld_list, drag_wetted_list, label="Drag based on wetted area")
plt.plot(ld_list, drag_frontal_list, label="Drag based on frontal area")
plt.plot(round(2 * (l_fus / (b_fus + h_fus)), 2), fuselage_drag(l_fus, b_fus, h_fus, rho_cruise), 'o', label="Actual fineness-ratio")
plt.axvline(x=ld_list[np.argmin(drag_wetted_list)], color='r', linestyle='--', label="Optimal fineness-ratio")
plt.xlabel("$l/d$ [-]")
plt.ylabel("Drag [N]")
plt.legend()
plt.xlim(0, 10)
plt.ylim(80, 250)

# plt.figure(4)
# plt.plot(ld_list, drag_perp_list)
# plt.axvline(x=ld_list[np.argmin(drag_perp_list)], color='r', linestyle='--', label="Optimal fineness-ratio")
# plt.xlabel("Fineness-ratio l/d [-]")
# plt.ylabel("Drag [N]")
# plt.legend()


difference_percentage_list = []
for i in range(len(drag_wetted_list)):
    difference_percentage = round(100 * (drag_wetted_list[i] - drag_frontal_list[i])/drag_wetted_list[i], 12)
    difference_percentage_list.append(difference_percentage)

print("Optimal fineness-ratio = ", str(round(ld_list[np.argmin(drag_wetted_list)], 1)))
print("Current fineness-ratio = ", str(round(2 * (l_fus / (b_fus + h_fus)),2)))
print("Maximum differerence drag estimation methods =", str(round(max(difference_percentage_list), 3)), "%")
print()

# plt.figure(2)
# plt.plot(ld_list, difference_percentage_list)


Re_list = np.arange(2.0*10**-1, 1*10**7, 1*10**3)
Cd_list = []

for i in range(len(Re_list)):
    Cd_list.append(Cd_Re(Re_list[i]))

plt.figure(5)
plt.semilogx(Re_list, Cd_list, label="Drag coefficient model")
plt.semilogx(Reynolds(d_bars, V_cruise/1.2), Cd_Re(Reynolds(d_bars, V_cruise/1.2)), "o", color="black", label=f"Cyl. with d = {d_bars} m and V = {round(V_cruise/1.2,2)} m/s")
plt.semilogx(Reynolds(d_bars, V_cruise), Cd_Re(Reynolds(d_bars, V_cruise)), "o", color="red", label=f"Cyl. with d = {d_bars} m and V = {round(V_cruise,2)} m/s")
plt.semilogx(Reynolds(d_hub, V_cruise/1.2), Cd_Re(Reynolds(d_hub, V_cruise/1.2)), "o", color="purple", label=f"Hub with d = {d_hub} m and V = {round(V_cruise/1.2,2)} m/s")
plt.semilogx(Reynolds(d_hub, V_cruise), Cd_Re(Reynolds(d_hub, V_cruise)), "o", color="orange", label=f"Hub with d = {d_hub} m and V = {round(V_cruise,2)} m/s")
plt.xlim(0.3*10**3, 1*10**7)
plt.ylim(0, 2)
plt.xlabel("Re [-]")
plt.ylabel("$C_D$ [-]")
plt.legend()

drag_vertical_flight = 0
drag_upper_structure_vertFlight = 0
for i in range(len(prop_structure_lenghts)):
    drag_upper_structure_vertFlight += cylinder_drag_perp(prop_structure_lenghts[i], d_bars, rho_cruise, 6)*4
drag_upper_structure_vertFlight += 2*cylinder_drag_perp(2.490, d_bars, rho_cruise, 6)
drag_upper_structure_vertFlight += 12*cylinder_drag_par(l_hub, d_hub, rho_cruise, 6)
drag_fuselage_vertFlight = Cd_Re(Reynolds(h_fus, 6)) * 0.5 * rho_cruise * 6 ** 2 * b_fus * l_fus
drag_skid_vert_bar_ver = cylinder_drag_par(0.5, d_bars, rho_cruise, V_cruise)
drag_skid_hor_bar_ver = cylinder_drag_perp(1.4, d_bars, rho_cruise, V_cruise)
drag_skid_vertFlight = 2*(drag_skid_vert_bar_ver*2 + drag_skid_hor_bar_ver)

drag_vertical_flight = drag_upper_structure_vertFlight + drag_fuselage_vertFlight + drag_skid_vertFlight
print("Total drag vertical flight = ", str(round(drag_vertical_flight, 3)), "N")

fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
drag_components = ["Upper structure", "Fuselage", "Skid"]
drag_components_values = [round(drag_upper_structure_vertFlight,3), round(drag_fuselage_vertFlight,3), round(drag_skid_vertFlight,3)]

wedges, texts, autotexts = ax.pie(drag_components_values, shadow=True, startangle=90, autopct=lambda pct: func(pct, drag_components_values),
                                  textprops=dict(color="w"))
ax.legend(wedges, drag_components,
          title="Drag Components",
          loc="center left",
          bbox_to_anchor=(0.9, 0, 0.5, 1))
plt.setp(autotexts, size=8, weight="bold")
ax.set_title(f"Total drag vertical flight = {round(drag_vertical_flight, 2)} N")

V_list = np.arange(0.1, 35, 0.2)
lower_indFactor = 0.8
higher_indFactor = 1.4
chosen_indFactor = 1.2
tot_drag_VSens_middle = []
tot_drag_VSens_lower = []
tot_drag_VSens_higher = []
tot_drag_VSens_chosen = []
for i in range(len(V_list)):
    tot_drag_VSens_middle.append(total_drag_calc(rho_cruise, V_list[i], 1))
    tot_drag_VSens_lower.append(total_drag_calc(rho_cruise, V_list[i], lower_indFactor))
    tot_drag_VSens_higher.append(total_drag_calc(rho_cruise, V_list[i], higher_indFactor))
    tot_drag_VSens_chosen.append(total_drag_calc(rho_cruise, V_list[i], chosen_indFactor))

plt.figure(7)
plt.fill_between(V_list*3.6, tot_drag_VSens_lower, tot_drag_VSens_higher, alpha=0.2)
plt.plot(V_list*3.6, tot_drag_VSens_higher, c='black', linestyle='--', label=f'$V$ higher = {higher_indFactor}$V_c$')
plt.plot(V_list*3.6, tot_drag_VSens_chosen, c='blue', linestyle='-.', label=f'$V$ chosen = {chosen_indFactor}$V_c$')
plt.plot(V_list*3.6, tot_drag_VSens_middle, label=f"$V_c$ = 100 km/h")
plt.plot(V_list*3.6, tot_drag_VSens_lower, c='black', linestyle='--', label=f'$V$ lower = {lower_indFactor}$V_c$')
plt.axvline(x=100, color='r', linestyle='--', label="Cruise velocity")
plt.axvline(x=(5*3.6), color='g', linestyle='--', label="Climb velocity")
plt.xlabel("$V$ [km/h]")
plt.ylabel("$D$ [N]")
plt.ylim(0, 500)
plt.legend()

d_test_verf = d_bars
drag_withSkinFr_list = []
drag_withoutSkinFr_list = []

for i in range(len(ld_list)):
    l_test_verf = ld_list[i] * d_test_verf
    if ld_list[i] < 2:
        C_D = ((1.17-0.81)/(0-2))*ld_list[i] + 1.17
    elif ld_list[i] > 2:
        C_D = 0.81

    Re = Reynolds(d_test_verf, V_cruise)
    friction_drag = (1.328 * 0.5*V_cruise**2 * l_test_verf * np.pi * d_test_verf)/(np.sqrt(Re))
    S0 = np.pi*(d_test_verf/2)**2
    drag_withSkinFr_list.append(C_D*0.5*rho_cruise*V_cruise**2 * S0 + friction_drag)
    drag_withoutSkinFr_list.append(C_D*0.5*rho_cruise*V_cruise**2 * S0)

plt.figure(8, figsize=(6, 2))
plt.plot(ld_list, drag_withSkinFr_list, label="Drag estimate including friction")
plt.plot(ld_list, drag_withoutSkinFr_list, label="Drag estimate without friction")
plt.legend()
plt.xlabel("$l/d$ [-]")
plt.ylabel("Drag [N]")

plt.show()

