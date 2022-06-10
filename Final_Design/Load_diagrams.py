import numpy as np
import matplotlib.pyplot as plt

class LoadDiagrams():

    def __init__(self, R, W, omega, P_a):
        self.rho = 1.225
        self.R_rotor = R
        self.n_rotors = 12
        self.W = W*9.81

        self.A_rotor = self.R_rotor**2 * np.pi
        self.A_T = self.A_rotor * 12
        self.disk_loading = self.W/self.A_T

        self.omega = omega*0.1047198
        self.v_h = np.sqrt(self.disk_loading/(2*self.rho))
        self.C_T = self.disk_loading/(self.rho*self.omega**2*self.R_rotor**2)
        self.V_C = 100/3.6
        self.P_h = 2*self.rho*self.A_T*self.v_h**3
        self.k = 1.15
        self.P_a = P_a
        self.FM = 0.7

    def get_drag(self, V):
        D_c = 350
        F = D_c / (0.5 * self.rho * self.V_C ** 2)
        D = F * 0.5 * self.rho * V ** 2
        return D

    def get_drag_horizontal(self, V_h):
        D_h = 60
        F = D_h / (0.5 * self.rho * 6 ** 2)
        D = F * 0.5 * self.rho * V_h ** 2
        return D

    def get_aoa(self, V):
        aoa = np.tanh(self.get_drag(V)/self.W)
        return aoa

    def get_parasitic_power(self, V):
        return self.get_drag(V)*V / self.FM

    def plot_aoa(self):
        V = np.arange(0, 100, 0.1)
        aoa = [self.get_aoa(i)*180/np.pi for i in V]
        plt.plot(V, aoa)
        plt.show()

    def get_advance_ratio(self, V, aoa = True):
        if aoa is True:
            AoA = self.get_aoa(V)
        else:
            AoA = aoa

        #Advance ratio in x-direction
        mu_x = V*np.cos(AoA)/(self.omega*self.R_rotor)

        return mu_x

    def get_inflow_ratio(self, V, aoa=True):
        if aoa is True:
            AoA = self.get_aoa(V)
        else:
            AoA = aoa

       #Get advance ratio in x and y direction
        mu_x = self.get_advance_ratio(V)

        #Set start inflow ratio to hover inflow ratio, then predict first inflow ration
        start_ratio = np.sqrt(self.C_T/2)
        inflow_ratio = mu_x * np.tan(AoA) + self.C_T / (2*np.sqrt(mu_x**2 + start_ratio**2))
       # list for keeping track of iterations
        while abs((inflow_ratio - start_ratio)/inflow_ratio) > 0.0005:
            start_ratio = inflow_ratio
            inflow_ratio =  mu_x * np.tan(AoA) + self.C_T / (2*np.sqrt(mu_x**2 + start_ratio**2))
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
        V = np.arange(0, 100/3.6, 0.1)
        T = [self.get_thrust(i, 0)/12 for i in V]
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

    def get_profile_power(self, V):
        mu = self.get_advance_ratio(V)
        K = 4.65
        Cd0 = 0.0085
        sigma = 0.1
        cp = (1+K*mu**2)*sigma*Cd0/8
        P = cp*self.A_T* (self.omega*self.R_rotor) ** 3
        return P / self.FM

    def get_power(self, V):
        P_i = self.get_v_inflow(V)*self.get_thrust(V, 0) / self.FM
        P_p = self.get_parasitic_power(V) / self.FM
        P_profile = self.get_profile_power(V)
        return P_i + P_p + P_profile


    def plot_power_actual(self):
        V = np.arange(0, 200/3.6, 0.01)
        P = [self.get_power(i)/12000 for i in V]
        P_i = [self.get_v_induced(i)*self.get_thrust(i, 0)/12000 for i in V]
        P_p = [self.get_advance_ratio(i)*self.omega*self.R_rotor*np.tan(self.get_aoa(i))*self.get_thrust(i, 0)/12000 for i in V]
        P_d = [self.get_parasitic_power(i)/12000 for i in V]
        P_prof = [self.get_profile_power(i)/12000 for i in V]
        P_a = np.ones(len(V))*self.P_a

        plt.plot(V*3.6, P_a, label='Available power', c='black', linestyle='--')
        plt.plot(V*3.6, P, label='Total power required ', c='black')
        plt.plot(V*3.6, P_i, label='Induced Power')
        plt.plot(V*3.6, P_p, label='Propulsive Power')
        plt.plot(V*3.6, P_d, label='Parasitic Power')
        plt.plot(V*3.6, P_prof, label='Profile Power')
        plt.xlabel('$V_\infty (km/h)$')
        plt.ylabel('P (kW)')
        plt.legend()
        plt.show()


    def performance(self):
        V = np.arange(0, 200/3.6, 0.1)
        P = np.array([self.get_power(i)/12000 for i in V])
        P_a = np.ones(len(V)) * self.P_a

        P_end = np.arange(0, min(P), 0.1)
        V_end = np.ones(len(P_end))*V[np.where(P==min(P))]

        slope_list = P/V
        slope = min(P/V)
        P_ran = np.arange(0, P[np.where(slope_list==slope)], 0.1)
        V_ran = np.ones(len(P_ran))*V[np.where(slope_list==slope)]
        range = slope*V

        index = 0
        for i, value in enumerate(P):
            if value > P_a[0]:
                index = i
                break

        P_max = np.arange(0, P_a[index], 0.1)
        V_max = np.ones(len(P_max))*V[index]

        plt.plot(V * 3.6, P_a, label='Available Power', c='black', linestyle='--')
        plt.plot(V * 3.6, P, label='Total Power Required', c='black')

        plt.plot(V_end*3.6, P_end, label=f'V max. endurance = {round(V_end[0]*3.6)} (km/h)')
        plt.plot(V_ran*3.6, P_ran, label=f'V max. range = {round(V_ran[0]*3.6)} (km/h)', c='orange')
        plt.plot(V_max*3.6, P_max, label=f'max. V attainable = {round(V_max[0]*3.6)} (km/h)', c='green')
        plt.plot(V * 3.6, range, linestyle=':', c='orange')

        plt.legend()
        plt.show()


    def plot_climb_performance(self):
        V = np.arange(0, 200/3.6, 0.1)
        V_c = [(15000*12-self.get_power(i)*self.k)/self.W for i in V]
        plt.plot(V*3.6, V_c)
        plt.show()


    def get_power_climb(self, V_c):
        P_c = ((V_c/(self.v_h*2)) + np.sqrt((V_c/(self.v_h*2))**2+1)) * self.P_h
        P_p = self.get_profile_power(0)
        P_d = V_c * self.get_drag_horizontal(V_c)
        return P_c + P_p + P_d

    def get_power_descent(self, V_c):
        P_d = (-0.125* V_c + 0.974 - 1.372*V_c**2 -1.718*V_c**3 - 0.655*V_c**4)*self.P_h
        P_p = self.get_profile_power(0)
        P_drag = V_c * self.get_drag_horizontal(V_c)
        return P_d + P_p + P_drag

    def get_thrust_climb(self, V_c):
        v_i = (-(V_c/(self.v_h*2)) + np.sqrt(((V_c/(self.v_h*2))**2)+1))*self.v_h
        T = 2*self.rho*self.A_T*(V_c+v_i)*v_i
        D = self.get_drag_horizontal(V_c)
        print(D)
        return T + D

    def get_inflow_climb(self, V_c):
        return (-(V_c/(self.v_h*2)) + np.sqrt(((V_c/(self.v_h*2))**2)+1))*self.v_h + V_c

    def plot_power_climb(self):
        V_c = np.arange(0, 10, 0.1)
        P = [self.get_power_climb(i) / 12000 for i in V_c]
        plt.plot(V_c, P)
        plt.xlabel('$V_{C}$ (m/s)')
        plt.ylabel('P (kW)')
        plt.show()

    def plot_thrust_climb(self):
        V_c = np.arange(0, 50, 0.1)
        T = [self.get_thrust_climb(i) for i in V_c]
        plt.plot(V_c, T)
        plt.show()

    def plot_climb_descent(self):
        V_c = np.arange(0, 10, 0.1)
        V_d = np.arange(-20, 0, 0.1) / self.v_h
        V_t = np.arange(-20, 20, 0.1)

        P_c = [self.get_power_climb(i) / 12000 for i in V_c]
        P_d = [self.get_power_descent(i)  / 12000 for i in V_d]
        P_a = np.ones(len(V_t))*self.P_a

        plt.plot(V_c, P_c, label='Normal working state')
        plt.plot(V_d*self.v_h, P_d, label='Vortex ring state')
        plt.plot(V_t, P_a, label='Power available', linestyle='--', c='black')
        plt.xlabel('$V_{c}$ and $V_{d}$')
        plt.ylabel('P (kW)')
        plt.legend()
        plt.show()


    def loads_cruise(self):
        V = np.arange(0, 100, 0.1)
        n = [self.get_thrust(i, 0)*np.cos(self.get_aoa(i)) / (self.W) for i in V]

        plt.plot(V, n)
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

    def program_input(self, V, V_c):
        print('For input:')
        print(f'Disk Loading = {self.disk_loading} (N/m^2)')
        print(f'Weight = {self.W} N')
        print(f'rpm = {self.omega/0.1047198}')
        print(f'AoA = {self.get_aoa(V)*180/np.pi} deg\n')
        print(f'The inputs to the program are:')
        print(f'Rotor diameter = {self.R_rotor*2} (m)')
        print(f'Inflow velocity = {self.get_v_inflow(V)} (m/s)')
        print(f'Thrust single engine = {self.get_thrust(V, 0)/12} (N)')
        print(f'Power single engine = {self.get_power(V)/12}(W)\n')
        print(f'Climb speed = {V_c} (m/s)')
        print(f'Thrust climb = {self.get_thrust_climb(V_c)/12} (N)')
        print(f'Power climb = {self.get_power_climb(V_c)/12} (N)')
        print(f'Inflow speed climb = {self.get_inflow_climb(V_c)} (m/s)')





# Input: Radius (m), Weight (kg), rpm, P_a
diagram = LoadDiagrams(0.9, 1256.8, 2011.5, 26.163)

diagram.program_input(0, 6)

#print(diagram.get_power_climb(0))
#print(diagram.get_power(0))
diagram.plot_power_actual()
#diagram.plot_climb_descent()
#diagram.plot_climb_performance()
diagram.performance()







