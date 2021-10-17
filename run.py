import numpy as np
import time
from matplotlib import pyplot as plt
from numpy.core.fromnumeric import var
import multiprocessing as mp
import csv

from tank_temp import Tank, Lamp, Thermals
from n2o_class import N2O

# tank data
k_carbon = 6           # W/m/K
k_alu = 237             # W/m/K
t_carbon = 2.5e-3         # m
t_alu = 1.5e-3            # m
L = 4.7                 # m
d = 0.2778              # m
S = np.pi*d*L           # m^2             
p_max = 6000            # kpa


# lamp data
n_lamp = 8
A_lamp = 4e-2 * 73e-2      # m^2
shape_factor = 0.00281 #0.00281  0.0029615889
emissivity = 0.6
lamp_temp = 899.39         # K

# sim data 
T = 2.5*3600#(57.53-5.09)*60                # s
dt = 0.001              # s  0.001 for stability if evap is considered 
P0 = 3750                # kPa



if __name__ == '__main__':

    tank = Tank(k_carbon, t_carbon, k_alu, t_alu, S)
    lamp = Lamp(emissivity, lamp_temp, shape_factor, n_lamp, A_lamp)
    n2o = N2O()     

    p_range = np.arange(3500, 5500+200, 200)
    times = np.ndarray((len(p_range),int(T/dt)))
    tank_pressures = np.ndarray((len(p_range),int(T/dt)))
    tank_temp = np.ndarray((len(p_range),int(T/dt)))


    for i in range(len(p_range)):

        T0 = n2o.get_temp(p_range[i])

        t0 = time.time()
        thermals = Thermals(tank, lamp, n2o, T, dt, T0, evap=True)
        thermals.forward_euler(p_max)
        t1 = time.time()

        print('run time: ', t1-t0, ' seconds') 

        times[i,:] = thermals.time
        tank_pressures[i,:] = thermals.P_tank
        tank_temp[i,:] = thermals.T_tank 

    
    with open('pressure.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["time","P_tank1","P_tank2","P_tank3","P_tank4","P_tank5","P_tank6","P_tank7","P_tank8","P_tank9","P_tank10"])
        for idx in range(int(T)):
            i = int(idx/dt)
            writer.writerow((times[0,i], tank_pressures[0,i],tank_pressures[1,i],tank_pressures[2,i],tank_pressures[3,i],tank_pressures[4,i],tank_pressures[5,i],tank_pressures[6,i],tank_pressures[7,i],tank_pressures[8,i],tank_pressures[9,i]))


    #var_pressure = np.array([37.5,37.9,39,39.75,40.7,41.07,41.5,41.9,42.53,42.8,43.1,43.6,44.1,44.26,44.6,45.47,45.53])*100
    #var_time = (np.array([5.09,7,13.27,17.46,23.24,25.49,28.42,30.59,35.17,37.22,38.58,42.35,46.22,47.53,50.21,57.25,57.53])-5.09)*60
    '''
    for i in range(len(p_range)):

        plt.plot(times[i], tank_pressures[i], label='P0 = '+str(p_range[i])+' kPa')
    #plt.scatter(var_time, var_pressure)

    plt.xlabel('time [s]')
    plt.ylabel('pressure [kPa]')
    plt.title('SIV Heating Time with 8 lamps')
    plt.legend(loc='best')
    plt.show()
    plt.savefig('siv_heating_time_parametric.png')
    '''
    