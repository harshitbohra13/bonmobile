import numpy as np
import matplotlib.pyplot as plt

class LoadDiagrams():

    def __init__(self, disk_loading, W, omega, AoA):
        self.AoA = AoA*np.pi/180
        self.rho = 1.225
        self.disk_loading = disk_loading
        self.n_rotors = 12
        self.W = W*9.81

        self.A_T = self.W / self.disk_loading
        self.A_rotor = self.A_T/self.n_rotors
        self.R_rotor = np.sqrt(self.A_rotor/np.pi)

        self.omega = omega*0.1047198
        self.v_h = np.sqrt(self.disk_loading/(2*self.rho))
        self.C_T = self.disk_loading/(self.rho*self.omega**2*self.R_rotor**2)
        self.V_C = 100/3.6

        self.Cdo = 0.86
        self.S = 100

    def print_rotor_data(self):
        print(f'Area rotor = {self.A_rotor}')
        print(f'Rotor diameter= {2*self.R_rotor}')
        print(f'Hover inflow speed= {self.v_h}')
        print(f'Thrust coefficient= {self.C_T}')

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

    def plot_inflow_speed(self):
        hover_inflow = np.sqrt(self.C_T/2)
        speed_range = np.arange(0, 50, 1)
        mu_x_lambda_h = [self.get_advance_ratio(i, self.AoA)[0]/hover_inflow for i in speed_range]
        print(mu_x_lambda_h)
        lambda_lambda_h = [self.get_inflow_ratio(i, self.AoA)/hover_inflow for i in speed_range]
        print(lambda_lambda_h)
        plt.plot(mu_x_lambda_h, lambda_lambda_h)
        plt.show()

    def get_v_i(self, V_c, AoA):
        v_i = (self.get_inflow_ratio(V_c, AoA) * self.omega * self.R_rotor) - (V_c * np.sin(AoA))
        return v_i

    def get_thrust(self, V_c, V_v):
        v_i = self.get_v_i(V_c, self.AoA)
        T = 2*self.rho*self.A_T * (v_i + V_v) * np.sqrt(V_c**2 + (V_c*np.sin(self.AoA) +v_i)**2)
        return T

    def get_load_factor(self, V_c, V_v):
        return self.get_thrust(V_c, V_v) / self.W

    def loads_cruise(self):
        x = np.arange(0, 40, 0.1)
        y = [self.get_load_factor(V_C, 0) for V_C in x]
        plt.plot(x, y)
        plt.xlabel('$V_{\inf} (m/s)$')
        plt.ylabel('n')
        plt.legend()
        plt.show()

    def loads_verticlimb(self):
        x = np.arange(0, 10, 0.1)
        y = [self.get_load_factor(0, V_v) for V_v in x]
        plt.plot(x, y)
        plt.xlabel('$V_{v} (m/s)$')
        plt.ylabel('n')
        plt.legend()
        plt.show()

    def loads_vertidescent(self):
        x = np.arange(0, -10, 0.1)
        y = [self.get_load_factor(0, V_v) for V_v in x]
        plt.plot(x, y)
        plt.xlabel('$V_{v} (m/s)$')
        plt.ylabel('n')
        plt.legend()
        plt.show()

    def inflow_thrust(self):
        v_c = np.arange(0, 36.4, 0.1)
        v_i = [self.get_v_i(V_h, self.AoA) for V_h in v_c]
        inflow_speed = v_i + v_c*np.sin(self.AoA)
        T = [self.get_thrust(V_h, 0) for V_h in v_c]
        plt.plot(inflow_speed, T)
        plt.xlabel('$V_{inflow} (m/s)$')
        plt.ylabel('T')
        plt.legend()
        plt.show()

    def inflow_aoa(self):
        AoA_range = np.arange(0, 15, 0.1)*np.pi/180
        v_i = [self.get_v_i(self.V_C, AoA)for AoA in AoA_range]
        plt.plot(np.arange(0, 15, 0.1), v_i)
        plt.xlabel('\u03B1 (\u03B1)')
        plt.ylabel('$V_{inflow} (m/s)$')
        plt.legend()
        plt.show()

    def inflow_cruise(self):
        V_c_range = np.arange(0, 200, 0.1)
        v_i = [self.get_v_i(V_c, self.AoA)for V_c in V_c_range]
        plt.plot(V_c_range, v_i)
        plt.xlabel('$V_{c} (m/s)$')
        plt.ylabel('$V_{inflow} (m/s)$')
        plt.legend()
        plt.show()

    def program_input(self):
        print('For input:')
        print(f'Disk Loading = {self.disk_loading} (N/m^2)')
        print(f'Weight = {self.W} kg')
        print(f'rpm = {self.omega/0.1047198}')
        print(f'AoA = {self.AoA*180/np.pi} deg\n')
        print(f'The inputs to the program are:')
        print(f'Rotor diameter = {self.R_rotor*2} (m)')
        print(f'Inflow velocity = {self.get_v_i(self.V_C, self.AoA)} (m/s)')
        print(f'Thrust single engine = {self.get_thrust(self.V_C, 0)/12} (N)')




# Input: Disk Loading (N/m^2), Weight (kg), rpm, AoA (deg)
diagram = LoadDiagrams(500, 1000, 6500, 0)
diagram.program_input()

