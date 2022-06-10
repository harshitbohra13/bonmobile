import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class PowerSensitivity():

    def get_power(self, V, R, W, omega, rho, D_c):
        n_rotors = 12
        W = 9.81*W
        R_rotor = R

        A_rotor = R_rotor ** 2 * np.pi
        A_T = A_rotor * 12
        disk_loading = W / A_T


        omega = omega * 0.1047198
        v_h = np.sqrt(disk_loading / (2 * rho))
        C_T = disk_loading / (rho * omega ** 2 * R_rotor ** 2)
        P_h = 2 * rho * A_T * v_h ** 3
        F = D_c / (0.5 * rho * (100/3.6) ** 2)
        D = F* 0.5* rho * V**2
        aoa = np.tanh(D / W)
        k = 1.15
        K_profile = 4.65
        Cd0 = 0.0085
        sigma = 0.1


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
        P_induced = v_inflow*T
        P_parasitic = D*V
        cp_profile = (1 + K_profile * mu_x ** 2) * sigma * Cd0 / 8
        P_profile = cp_profile * A_T * omega ** 3 * R_rotor ** 3
        P = P_induced + P_parasitic + P_profile

        return P, T

    def power_thrust(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 100, 0.1)
        P_T = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/self.get_power(i, disk_loading, W, omega, rho, D_c)[0] for i in V]
        plt.plot(V, P_T)
        plt.show()

    def disk_loading_P(self,disk_loading, W, omega, rho, D_c, P_a):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        P_actual = np.ones(len(V)) * P_a

        lower = [self.get_power(i, disk_loading*0.5, W, omega, rho, D_c)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading*1.5, W, omega, rho, D_c)[0]/12000 for i in V]

        plt.plot(V * 3.6, P_actual, c='r', label='Available power')
        plt.plot(V * 3.6, lower, c='black', linestyle='--', label=f'Disk loading = {disk_loading*0.5}')
        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'Disk loading = {disk_loading*1.5}')
        plt.plot(V*3.6, P, c='b', label=f'Disk loading = {disk_loading}')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def power_radius(self, R, W, omega, rho, D_c, P_a):
        rpms = np.arange(1000, 4000, 500)
        R = np.arange(0.5*R, 1.5*R, 0.01)
        R_optimal_list = []

        for rpm in rpms:
             P = [self.get_power(0, i, W, rpm, rho, D_c)[0]/(12000*0.6) for i in R]
             min_P = min(P)
             R_optimal = R[P.index(min_P)]
             R_optimal_list.append(R_optimal)
             plt.plot(R, P, label=f'rpm={rpm}')
             plt.plot(R_optimal, min_P, marker='o', c='r')

        P_actual = np.ones(len(R)) * P_a
        plt.plot(R, P_actual, c='black', linestyle='--', label='Available power')
        plt.xlabel('Rotor radius')
        plt.ylim([15,60])
        plt.ylabel('Power required cruise')
        plt.legend()
        plt.show()

        plt.plot(rpms, R_optimal_list)
        plt.xlabel('rpm')
        plt.ylabel('R rotor')
        plt.show()

    def power_rpm_final(self, R, W, omega, rho, D_c, P_a):
        rpms = np.arange(0, 4000, 100)
        radius = np.arange(0.7, 1.1, 0.1)
        rpm_optimal_list = []

        for R in radius:
            P = [self.get_power(0, R, W, i, rho, D_c)[0] / (12000 * 0.7) for i in rpms]
            min_P = min(P)
            rpm_optimal = rpms[P.index(min_P)]
            rpm_optimal_list.append(rpm_optimal)
            plt.plot(rpms, P, label=f'R={round(R, 1)}')
            plt.plot(rpm_optimal, min_P, marker='o', c='r')

        P_actual = np.ones(len(rpms)) * P_a
        #plt.plot(rpms, P_actual, c='black', linestyle='--', label='Available power')
        plt.xlabel('rpm')
        plt.ylim([10, 60])
        plt.ylabel('Power required cruise')

        bat1 = (22/5000)*rpms
        bat2 = (32/4000)*rpms
        bat3 = (52/4000)*rpms

        plt.plot(rpms, bat1, label='bat 1', linestyle =':')
        plt.plot(rpms, bat2, label='bat 2', linestyle =':')
        plt.plot(rpms, bat3, label='bat 3', linestyle =':')
        plt.legend()
        plt.show()

    def rpm_velocity(self, R, W, omega, rho, D_c, P_a):
        V_c = np.arange(0, 140, 0.1)

    def power_rpm_hover(self, R, W, omega, rho, D_c, P_a):
        rpms = np.arange(0, 4100, 100)
        P = np.array([self.get_power(0, R, W, i, rho, D_c)[0] / (12000*0.7) for i in rpms])
        bat = (52 / 4000) * rpms

        P_a = bat - P

        plt.plot(rpms, P, label='Power required')
        plt.plot(rpms, P_a, label='Excess power')
        plt.plot(rpms, bat, label='Power available')
        plt.ylim([0, 50])
        plt.xlim([0, 4000])
        plt.xlabel('\u03C9 (rpm)')
        plt.ylabel('P (kW)')
        plt.legend()
        plt.show()


    def D_c_P(self, disk_loading, W, omega, rho, D_c, P_a):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        P_actual = np.ones(len(V))* P_a

        lower = [self.get_power(i, disk_loading, W, omega, rho, D_c*0.5)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c*1.5)[0]/12000 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, P_actual, c='r', label='Available power')
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'D lower = {D_c*0.5} (N)')
        plt.plot(V*3.6, P, c='b', label=f'D = {D_c} (N)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D upper = {D_c*1.5} (N)')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def omega_P(self, disk_loading, W, omega, rho, D_c, P_a):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        P_actual = np.ones(len(V))* P_a

        lower = [self.get_power(i, disk_loading, W, omega*0.5, rho, D_c)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega*1.5, rho, D_c)[0]/12000 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, P_actual, c='r', label='Available power')
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'\u03A9 lower = {omega*0.5} (rpm)')
        plt.plot(V*3.6, P, c='b', label=f'D = {omega} (N)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'\u03A9 upper = {omega*1.5} (rpm)')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()



    def rho_P(self, disk_loading, W, omega, rho, D_c, P_a):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        P_actual = np.ones(len(V)) * P_a

        lower = [self.get_power(i, disk_loading, W, omega, 1.1731597739917108, D_c)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, P, c='b', label=f'\u03C1 = {rho} (kg/m^3)')
        plt.plot(V * 3.6, P_actual, c='r', label='Available power')
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'\u03C1 lower = {rho*0.5} (kg/m^3)')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def W_P(self, disk_loading, W, omega, rho, D_c, P_a):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[0]/12000 for i in V]
        P_actual = np.ones(len(V)) * P_a

        lower = [self.get_power(i, disk_loading, W*0.5, omega, rho, D_c)[0]/12000 for i in V]
        upper = [self.get_power(i, disk_loading, W*1.5, omega, rho, D_c)[0]/12000 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'W lower = {W*0.5} (kg)')
        plt.plot(V * 3.6, P_actual, c='r', label='Available power')
        plt.plot(V*3.6, P, c='b', label=f'D = {W} (kg)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D lower = {W*1.5} (kg)')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$P$ (kW)')
        plt.legend()
        plt.show()

    def disk_loading_T(self,disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]
        lower = [self.get_power(i, disk_loading*0.5, W, omega, rho, D_c)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading*1.5, W, omega, rho, D_c)[1]/12 for i in V]

        plt.plot(V * 3.6, lower, c='black', linestyle='--', label=f'Disk loading = {disk_loading*0.5}')
        plt.fill_between(V*3.6, lower, upper, alpha=0.2, label=f'Disk loading = {disk_loading}')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'Disk loading = {disk_loading*1.5}')
        plt.plot(V*3.6, P, c='r')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def D_c_T(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 200/3.6, 1)
        lower = [self.get_power(i, disk_loading, W, omega, rho, D_c*0.5)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c*1.5)[1]/12 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'D lower = {D_c*0.5} (N)')
        plt.plot(V*3.6, P, c='r', label=f'D = {D_c} (N)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D upper = {D_c*1.5} (N)')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def rho_T(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]
        lower = [self.get_power(i, disk_loading, W, omega, rho*0.5, D_c)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, P, c='r', label=f' \u03C1 = {rho}($kg/m^3$')
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'\u03C1 = {rho*0.5}(N)')
        plt.xlabel('$V_{\infty} (km/h)$')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def W_T(self, disk_loading, W, omega, rho, D_c):
        V = np.arange(0, 200/3.6, 1)

        P = [self.get_power(i, disk_loading, W, omega, rho, D_c)[1]/12 for i in V]
        lower = [self.get_power(i, disk_loading, W*0.5, omega, rho, D_c)[1]/12 for i in V]
        upper = [self.get_power(i, disk_loading, W*1.5, omega, rho, D_c)[1]/12 for i in V]

        plt.fill_between(V*3.6, lower, upper, alpha=0.2)
        plt.plot(V*3.6, lower, c='black', linestyle='--', label=f'W lower = {W*0.5} (kg)')
        plt.plot(V*3.6, P, c='r', label=f'D = {W} (kg)')
        plt.plot(V*3.6, upper, c='black', linestyle='--', label=f'D lower = {W*1.5} (kg)')
        plt.xlabel('$V_{\infty}$(km/h)')
        plt.ylabel('$T$ (N)')
        plt.legend()
        plt.show()

    def do_sensitivity_analysis(self, disk_loading, W, omega, rho, D_c, P_a):
        self.power_radius(disk_loading, W, omega, rho, D_c, P_a)
        self.power_rpm_final(disk_loading, W, omega, rho, D_c, P_a)
        self.power_rpm_hover(disk_loading, W, omega, rho, D_c, P_a)
        #self.D_c_P(disk_loading, W, omega, rho, D_c, P_a)
        #self.disk_loading_P(disk_loading, W, omega, rho, D_c, P_a)
        #self.rho_P(disk_loading, W, omega, rho, D_c, P_a)
        #self.W_P(disk_loading, W, omega, rho, D_c, P_a)
        #self.omega_P(disk_loading, W, omega, rho, D_c, P_a)

        #self.D_c_T(disk_loading, W, omega, rho, D_c)
        #self.disk_loading_T(disk_loading, W, omega, rho, D_c)
        #self.rho_T(disk_loading, W, omega, rho, D_c)
        #self.W_T(disk_loading, W, omega, rho, D_c)


#R, W, omega, rho, D_c
sens = PowerSensitivity()
sens.do_sensitivity_analysis(0.9, 1200, 2011.5, 1.225, 350, 29)
#sens.power_thrust(500, 1500, 5000, 1.225, 400)