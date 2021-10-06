import numpy as np
from matplotlib import pyplot as plt

class Tank():
    def __init__(self, k_carbon, t_carbon, k_alu, t_alu, A):
        self.k_c = k_carbon
        self.t_c = t_carbon
        self.k_a = k_alu
        self.t_a = t_alu
        self.A = A
        self.m = 174       # kg
        self.V = 0.259      # m3

class Lamp():
    def __init__(self, eps, Temp, shape_factor):
        self.eps = eps
        self.T = Temp
        self.shape_factor = shape_factor
        self.sig = 5.67e-8

# dummy class 
class N2O():
    def __init__(self):
        self.rho = 300
        self.cp = 1000


class Thermals():
    def __init__(self, tank, lamp, n2o, T, dt, T_amb, evap=False):
        self.tank = tank
        self.lamp = lamp
        self.n2o = n2o
        self.T = T
        self.dt = dt 
        self.T_amb = T_amb
        self.evap = evap

    def thermal_res(self, T):
        h_r = 4 * self.lamp.eps * self.lamp.sig * self.lamp.shape_factor * ((T + self.lamp.T)/2)**3 
        R = (1/h_r + self.tank.t_c/self.tank.k_c + self.tank.t_a/self.tank.k_a) * 1/self.tank.A
        
        return R

    def get_fluid_properties(self, T):
        vap_frac = self.n2o.vap_mass_frac(T)            # mass fraction of vapor for mass averaged density and Cp
        #mass_frac = self.n2o.m_V(T)/self.n2o.m_L(T)
        #print(vap_frac)
        rho = (1-vap_frac) * self.n2o.rho_V(T) + vap_frac * self.n2o.rho_L(T) 
        cp = (1-vap_frac) * self.n2o.cp_V(T) + vap_frac * self.n2o.cp_L(T) 
        
        return rho, cp, vap_frac


    def forward_euler(self, cutoff):
        self.time = np.arange(0, self.T, self.dt)
        self.T_tank = np.ones(len(self.time))*self.T_amb
        self.rho_n2o = np.ndarray(len(self.time))
        self.cp_n2o = np.ndarray(len(self.time))
        self.vap_frac = np.ndarray(len(self.time))

        self.T_tank[0] = self.T_amb
        self.rho_n2o[0], self.cp_n2o[0], self.vap_frac[0] = self.get_fluid_properties(self.T_tank[0])

        for i in range(1, len(self.time)):
            R = self.thermal_res(self.T_tank[i-1])
            deltaT = self.lamp.T - self.T_tank[i-1]

            self.rho_n2o[i], self.cp_n2o[i], self.vap_frac[i] = self.get_fluid_properties(self.T_tank[i-1])
            
            evap_mass = (self.vap_frac[i] - self.vap_frac[i-1]) * self.tank.m

            if self.evap:
                Q_evap = self.n2o.hvap * evap_mass 

            else:
                Q_evap = 0

            #Q_conv = self.conv_heat_transfer(self.T_tank[i-1])

            self.T_tank[i] = self.T_tank[i-1] + self.dt * ((deltaT/R  - Q_evap) * 1 / self.rho_n2o[i] / self.cp_n2o[i] / self.tank.V)

            if self.n2o.p(self.T_tank[i]) > cutoff:
                print('max pressure reached after: ', self.time[i], ' seconds')
                break

