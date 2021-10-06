import numpy as np
import time
from matplotlib import pyplot as plt

from tank_temp import Tank, Lamp, Thermals
from n2o_class import N2O

# tank data
k_carbon = 10           # W/m/K
k_alu = 200             # W/m/K
t_carbon = 2e-3         # m
t_alu = 3e-3            # m
A_tank = 4              # m^2
p_max = 7000            # kpa


# lamp data
shape_factor = 0.00269
emissivity = 0.6
lamp_temp = 350         # K

# sim data 
T = 3600                # s
dt = 0.001              # s  0.001 for stability if evap is considered 
T0 = 288                # K


if __name__ == '__main__':

    tank = Tank(k_carbon, t_carbon, k_alu, t_alu, A_tank)
    lamp = Lamp(emissivity, lamp_temp, shape_factor)
    n2o = N2O()             

    t0 = time.time()
    thermals = Thermals(tank, lamp, n2o, T, dt, T0, evap=True)
    thermals.forward_euler(p_max)
    t1 = time.time()

    print('run time: ', t1-t0, ' seconds')


    plt.plot(thermals.time, thermals.T_tank)
    plt.xlabel('time [s]')
    plt.ylabel('temp [K]')
    plt.show()