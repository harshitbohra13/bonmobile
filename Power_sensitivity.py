import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class PowerSensitivity():

    def get_power(self, V, disk_loading, W, omega, rho, D_c):
        n_rotors = 12
        W = 9.81*W

        A_T = W / disk_loading
        A_rotor = A_T / n_rotors
        R_rotor = np.sqrt(A_rotor / np.pi)
        omega = omega * 0.1047198
        v_h = np.sqrt(disk_loading / (2 * rho))
        C_T = disk_loading / (rho * omega ** 2 * R_rotor ** 2)
        P_h = 2 * rho * A_T * v_h ** 3
        F = D_c / (0.5 * rho * (100/3.6) ** 2)
        D = F* 0.5* rho * V**2
        aoa = np.tanh(D / W)

        mu_x = V * np.cos(aoa) / (omega * R_rotor)

        # Set start inflow ratio to hover inflow ratio, then predict first inflow ration
        start_ratio = np.sqrt(C_T / 2)
        inflow_ratio = mu_x * np.tan(aoa) + C_T / (2 * np.sqrt(mu_x ** 2 + start_ratio ** 2))

        # list for keeping track of iterations
        while abs((inflow_ratio - start_ratio) / inflow_ratio) > 0.0005:
            start_ratio = inflow_ratio
            inflow_ratio = mu_x * np.tan(aoa) + C_T / (2 * np.sqrt(mu_x ** 2 + start_ratio ** 2))

        v_inflow = inflow_ratio * omega * R_rotor
        v_induced = inflow_ratio * omega * R_rotor - V * np.sin(aoa)

        T = 2 * rho * A_T * v_induced * np.sqrt(V ** 2 + (v_inflow ** 2))
        P = v_inflow*T

        return P, T

    def power_thrust(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 100, 0.1)
        P_T = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/self.get_power(i, disk_loading, W, omega, rho, D_c)[0] for i in V]
        plt.plot(V, P_T)
        plt.show()


    def disk_loading_P(self,disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        lower = [self.get_power(i, disk_loading*0.5, W, omega, rho, D_c)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading*1.5, W, omega, rho, D_c)[0]/12000 for i in V]

        plt.plot(V * 3.6, lower, c='black', linestyle='--', label=f'Disk loading = {disk_loading*0.5}')
        plt.fill_between(V*3.6, lower, upper, alpha=0.2, label=f'Disk loading = {disk_loading}')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'Disk loading = {disk_loading*1.5}')
        plt.plot(V*3.6, P, c='r')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def D_c_P(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        lower = [self.get_power(i, disk_loading, W, omega, rho, D_c*0.5)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c*1.5)[0]/12000 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'D lower = {D_c*0.5} (N)')
        plt.plot(V*3.6, P, c='r', label=f'D = {D_c} (N)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D upper = {D_c*1.5} (N)')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def rho_P(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        lower = [self.get_power(i, disk_loading, W, omega, rho*0.5, D_c)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, P, c='r', label=f'\u03C1 = {rho} (kg/m^3)')
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'\u03C1 lower = {rho*0.5} (kg/m^3)')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def W_P(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        lower = [self.get_power(i, disk_loading, W*0.5, omega, rho, D_c)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading, W*1.5, omega, rho, D_c)[0]/12000 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'W lower = {W*0.5} (kg)')
        plt.plot(V*3.6, P, c='r', label=f'D = {W} (kg)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D lower = {W*1.5} (kg)')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def disk_loading_T(self,disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]
        lower = [self.get_power(i, disk_loading*0.5, W, omega, rho, D_c)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading*1.5, W, omega, rho, D_c)[1]/12 for i in V]

        plt.plot(V * 3.6, lower, c='black', linestyle='--', label=f'Disk loading = {disk_loading*0.5}')
        plt.fill_between(V*3.6, lower, upper, alpha=0.2, label=f'Disk loading = {disk_loading}')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'Disk loading = {disk_loading*1.5}')
        plt.plot(V*3.6, P, c='r')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def D_c_T(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]
        lower = [self.get_power(i, disk_loading, W, omega, rho, D_c*0.5)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c*1.5)[1]/12 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'D lower = {D_c*0.5} (N)')
        plt.plot(V*3.6, P, c='r', label=f'D = {D_c} (N)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D upper = {D_c*1.5} (N)')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def rho_T(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]
        lower = [self.get_power(i, disk_loading, W, omega, rho*0.5, D_c)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, P, c='r', label=f' \u03C1 = {rho}($kg/m^3$')
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'\u03C1 = {rho*0.5}(N)')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def W_T(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 140/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]
        lower = [self.get_power(i, disk_loading, W*0.5, omega, rho, D_c)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading, W*1.5, omega, rho, D_c)[1]/12 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'W lower = {W*0.5} (kg)')
        plt.plot(V*3.6, P, c='r', label=f'D = {W} (kg)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D lower = {W*1.5} (kg)')
        plt.xlabel('$V_{C}$(km/h)')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def do_sensitivity_analysis(self, disk_loading, W, omega, rho, D_c):
        self.D_c_P(disk_loading, W, omega, rho, D_c)
        self.disk_loading_P(disk_loading, W, omega, rho, D_c)
        self.rho_P(disk_loading, W, omega, rho, D_c)
        self.W_P(disk_loading, W, omega, rho, D_c)
        self.D_c_T(disk_loading, W, omega, rho, D_c)
        self.disk_loading_T(disk_loading, W, omega, rho, D_c)
        self.rho_T(disk_loading, W, omega, rho, D_c)
        self.W_T(disk_loading, W, omega, rho, D_c)


#disk_loading, W, omega, rho, D_c
sens = PowerSensitivity()
#sens.do_sensitivity_analysis(500, 1500, 2000, 1.225, 400)
sens.power_thrust(500, 1500, 5000, 1.225, 400)