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
        self.P_h = 2*self.rho*self.A_T*self.v_h**3

        self.Cdo = 0.86
        self.S = 100


    def print_rotor_data(self):
        print(f'Area rotor = {self.A_rotor}')
        print(f'Rotor diameter= {2*self.R_rotor}')
        print(f'Hover inflow speed= {self.v_h}')
        print(f'Thrust coefficient= {self.C_T}')

    def get_aoa(self, V):
        D_c = 350
        F = D_c/(0.5*self.rho*self.V_C**2)
        D  = F* 0.5*self.rho * V**2
        aoa = np.tanh(D/self.W)
        return aoa

    def get_advance_ratio(self, V, aoa = True):
        if aoa is True:
            AoA = self.get_aoa(V)
        else:
            AoA = aoa

        #Advance ratio in x-direction
        mu_x = V*np.cos(AoA)/(self.omega*self.R_rotor)

        #Advance ratio in y-direction
        mu_y = V*np.sin(AoA)/(self.omega*self.R_rotor)

        return mu_x, mu_y

    def get_inflow_ratio(self, V, aoa=True):
        if aoa is True:
            AoA = self.get_aoa(V)
        else:
            AoA = aoa

       #Get advance ratio in x and y direction
        mu_x, mu_y = self.get_advance_ratio(V)

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

        alpha_list = np.array([-2, 0, 2, 4, 6, 8, 10])*np.pi/180

        for aoa in alpha_list:
            mu_x_lambda_h = [self.get_advance_ratio(i, aoa)[0]/hover_inflow for i in speed_range]
            lambda_lambda_h = [self.get_inflow_ratio(i, aoa)/hover_inflow for i in speed_range]
            plt.plot(mu_x_lambda_h, lambda_lambda_h, label = f'\u03B1 = {aoa*180/np.pi}\u00b0')

        plt.xlabel('$\u03BB/\u03BB_{h}$')
        plt.ylabel('\u03BC/\u03BB_{h}')
        plt.legend()
        plt.show()

    def get_v_inflow(self, V, aoa=True):
        if aoa is True:
            AoA = self.get_aoa(V)
        else:
            AoA = aoa

        v_inflow = self.get_inflow_ratio(V, AoA) * self.omega * self.R_rotor
        return v_inflow

    def get_v_induced(self, V):
        v_induced= self.get_inflow_ratio(V) * self.omega * self.R_rotor - V*np.sin(self.get_aoa(V))
        return v_induced

    def plot_induced(self):
        v_c = np.arange(0, 100, 0.1)
        v_induced = [self.get_v_induced(V) for V in v_c]
        plt.plot(v_c, v_induced)
        plt.xlabel('V_{c} (m/s)')
        plt.ylabel('V_{i} (m/s)')
        plt.legend()
        plt.show()

    def get_thrust(self, V_h, V_v):
        T = 2*self.rho*self.A_T * (self.get_v_induced(V_h) + V_v) * np.sqrt(V_h**2 + (self.get_v_inflow(V_h))**2)
        return T

    def plot_thrust(self):
        V = np.arange(0, 40, 0.1)
        T = [self.get_thrust(i, 0) for i in V]
        plt.plot(V, T)
        plt.xlabel('$V_c {m/s}}$')
        plt.ylabel('T (N)')
        plt.legend()
        plt.show()

    def plot_power_curves(self):
        speed_range = np.arange(0, 36, 1)

        alpha_list = np.array([-2, 0, 2, 4, 6, 8, 10])*np.pi/180

        for aoa in alpha_list:
            P = [self.get_v_inflow(i, aoa)*self.get_thrust(i, 0)/12 for i in speed_range]
            plt.plot(speed_range, P, label = f'\u03B1 = {aoa*180/np.pi}\u00b0')

        plt.xlabel('$V_c {m/s}}$')
        plt.ylabel('P (W)')
        plt.legend()
        plt.show()

    def get_power(self, V):
        return self.get_v_inflow(V)*self.get_thrust(V, 0)

    def plot_power(self):
        V = np.arange(0, 40, 0.1)
        P_T = [self.get_v_inflow(i)*self.get_thrust(i, 0) for i in V]
        P_i = [self.get_v_induced(i)*self.get_thrust(i, 0) for i in V]
        P_p = [self.get_advance_ratio(i)*np.tan(self.get_aoa(i))*self.get_thrust(i, 0) for i in V]

        plt.plot(V, P_T)
        plt.plot(V, P_i)
        plt.plot(V, P_p)
        plt.xlabel('$V_c {m/s}}$')
        plt.ylabel('P (W)')
        plt.legend()
        plt.show()

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
        v_i = [self.get_v_inflow(V_h, self.AoA) for V_h in v_c]
        inflow_speed = v_i + v_c*np.sin(self.AoA)
        T = [self.get_thrust(V_h, 0) for V_h in v_c]
        plt.plot(inflow_speed, T)
        plt.xlabel('$V_{inflow} (m/s)$')
        plt.ylabel('T')
        plt.legend()
        plt.show()

    def inflow_aoa(self):
        AoA_range = np.arange(0, 15, 0.1)*np.pi/180
        v_i = [self.get_v_inflow(self.V_C, AoA)for AoA in AoA_range]
        plt.plot(np.arange(0, 15, 0.1), v_i)
        plt.xlabel('\u03B1 (\u03B1)')
        plt.ylabel('$V_{inflow} (m/s)$')
        plt.legend()
        plt.show()

    def inflow_cruise(self):
        V_c_range = np.arange(0, 36, 0.1)
        v_i = [self.get_v_inflow(V_c, self.AoA)for V_c in V_c_range]
        plt.plot(V_c_range, v_i)
        plt.xlabel('$V_{c} (m/s)$')
        plt.ylabel('$V_{inflow} (m/s)$')
        plt.legend()
        plt.show()

    def program_input(self):
        print('For input:')
        print(f'Disk Loading = {self.disk_loading} (N/m^2)')
        print(f'Weight = {self.W} N')
        print(f'rpm = {self.omega/0.1047198}')
        print(f'AoA = {self.get_aoa(0)*180/np.pi} deg\n')
        print(f'The inputs to the program are:')
        print(f'Rotor diameter = {self.R_rotor*2} (m)')
        print(f'Inflow velocity = {self.get_v_inflow(0)} (m/s)')
        print(f'Thrust single engine = {self.get_thrust(0, 0)/12} (N)')
        print(f'Power single engine = {self.get_power(0)/12}(W)')




# Input: Disk Loading (N/m^2), Weight (kg), rpm, AoA (deg)
diagram = LoadDiagrams(500, 1000, 2000, 8)

diagram.program_input()

diagram.plot_power()

