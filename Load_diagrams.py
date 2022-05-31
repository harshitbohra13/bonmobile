import numpy as np
import matplotlib.pyplot as plt

class LoadDiagrams():

    def __init__(self, disk_loading):
        self.rho = 1.225
        self.disk_loading = disk_loading
        self.n_rotors = 12
        self.R_rotor = 0.62
        self.A_rotor = self.R_rotor**2 * np.pi
        self.W = 4414
        self.A_T = self.A_rotor * self.n_rotors
        self.C_T = 0.0016385
        self.omega = 628.318531
        self.v_h = np.sqrt(self.disk_loading/2*self.rho)
        self.C_T =

    def get_advance_ratio(self, V_c, AoA):

        #Advance ratio in x-direction
        mu_x = V_c*np.cos(AoA)/(self.omega*self.R_rotor)

        #Advance ratio in y-direction
        mu_y = V_c*np.sin(AoA)/(self.omega*self.R_rotor)

        return mu_x, mu_y

    def get_inflow_ratio(self, V_c, AoA):
       #Get advance ratio in x and y direction
        mu_x, mu_y = self.get_advance_ratio(V_c, AoA)

        #Set start inflow ratio to hover inflow ratio, then predict first inflow ration
        start_ratio = np.sqrt(self.C_T/2)
        inflow_ratio = mu_x * np.tan(AoA) + self.C_T / (2*np.sqrt(mu_x**2 + start_ratio**2)) + mu_y
       # list for keeping track of iterations
        while abs((inflow_ratio - start_ratio)/inflow_ratio) > 0.0005:
            start_ratio = inflow_ratio
            inflow_ratio =  mu_x * np.tan(AoA) + self.C_T / (2*np.sqrt(mu_x**2 + start_ratio**2)) + mu_y
        return inflow_ratio

    def plot_inflow_speed(self, AoA):
        hover_inflow = np.sqrt(self.C_T/2)
        speed_range = np.arange(0, 50, 1)
        mu_x_lambda_h = [self.get_advance_ratio(i, AoA)[0]/hover_inflow for i in speed_range]
        print(mu_x_lambda_h)
        lambda_lambda_h = [self.get_inflow_ratio(i, AoA)/hover_inflow for i in speed_range]
        print(lambda_lambda_h)
        plt.plot(mu_x_lambda_h, lambda_lambda_h)
        plt.show()

    def get_trust(self, V_c, V_v, AoA):
        v_i = (self.get_inflow_ratio(V_c, AoA)*self.omega*self.R_rotor) - (V_c*np.sin(AoA))
        print(v_i)
        T = 2*self.rho*self.A_T * (v_i + V_v) * np.sqrt(V_c**2 + (V_c*np.sin(AoA) +v_i)**2)
        return T

    def get_load_factor(self, V_c, V_v, AoA):
        return self.get_trust(V_c, V_v, AoA) / self.W

    def get_load_diagram_cruise(self, AoA):
        x = np.arange(0, 38.9, 0.1)
        y = [self.get_load_factor(V_C, 0, AoA) for V_C in x]
        plt.plot(x, y)
        plt.xlabel('$V_{\inf} (m/s)$')
        plt.ylabel('n')
        plt.legend()
        plt.show()

    def get_load_diagram_ax_climb(self, AoA):
        x = np.arange(0, 10, 0.1)
        y = [self.get_load_factor(0, V_v, AoA) for V_v in x]
        plt.plot(x, y)
        plt.xlabel('$V_{v} (m/s)$')
        plt.ylabel('n')
        plt.legend()
        plt.show()

    def get_load_diagram_ax_descent(self, AoA):
        x = np.arange(0, -10, 0.1)
        y = [self.get_load_factor(0, V_v, AoA) for V_v in x]
        plt.plot(x, y)
        plt.xlabel('$V_{v} (m/s)$')
        plt.ylabel('n')
        plt.legend()
        plt.show()





diagram = LoadDiagrams(1.225)

diagram.get_load_diagram_cruise(3*np.pi/180)
diagram.get_load_diagram_ax_climb(3*np.pi/180)
#diagram.get_load_diagram_ax_descent(0.1)
#diagram.get_3d_load_diagram(0.1)
#diagram.plot_inflow_speed(0.08)
#print(diagram.get_advance_ratio(100, 0))
#print(diagram.get_inflow_ratio(100, 0))
print(diagram.get_trust(0, 0, 0))

# # If plot is needed plots the progress of iterations
# if plot is True:
#     print(iteration_list)
#     plt.plot(range(len(iteration_list)), iteration_list)
#     plt.show()

# iteration_list = [start_ratio]
#         if threeD is True:
#             inflow_ratio_list = []
#             for values in inflow_ratio[0]:
#                 inflow_ratio_value = values
#                 while abs((inflow_ratio_value - start_ratio)/inflow_ratio_value) > 0.0005:
#                     start_ratio = inflow_ratio_value
#                     inflow_ratio_value =  mu_x * np.tan(AoA) + self.C_T / (2*np.sqrt(mu_x**2 + start_ratio**2)) + mu_y
#             return inflow_ratio_list
#         #Iterate till error becomes acceptably small
#         else: