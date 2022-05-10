import numpy as np

#aluminum yield strenght 270 MPA

r_outer = 0.075 #[m]
r_inner = 0.070 #[m]
structure_length = 2 #[m]
g = 9.8
Force = 400*g #[N]
M = Force *structure_length
I = np.pi/4 * (r_outer**4-r_inner**4)#[moment of inertia]
max_stress = r_outer* M/(I)

print("max stress", max_stress/(10**6), "MPa")

#