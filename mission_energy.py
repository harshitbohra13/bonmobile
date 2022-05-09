from power_concept1 import P_hover, P_climb
#times

t_hover =
t_climb = 
t_descend =

t_cruise = 40/100



# Total energy
E =  2 * P_hover*t_hover + P_climb*t_climb + P_cruise*t_cruise + P_descend * t_descend
print("Total mission energy:", E, "J")