import scipy.optimize as opt

from run import k_carbon, k_alu, t_carbon, t_alu, S, p_max, emissivity, lamp_temp, n_lamp, A_lamp
from tank_temp import Tank, Lamp, Thermals
from n2o_class import N2O

P0 = 3750               # kPa
time = 60*(57.53-5.09)         # s
p_measured = 4553       # kPa
dt = 10                 # s
SF = 0.002734757460420952              # intial guess for shape factor
bounds = [0.000001, 0.005]

# find shape fator that corresponds with these measured values 

tank = Tank(k_carbon, t_carbon, k_alu, t_alu, S)
lamp = Lamp(emissivity, lamp_temp, SF, n_lamp, A_lamp)
n2o = N2O()    
T0 = n2o.get_temp(P0)  
print(T0)

def func(shape_factor):
    lamp = Lamp(emissivity, lamp_temp, shape_factor, n_lamp, A_lamp)
    thermals = Thermals(tank, lamp, n2o, time, dt, T0, evap=False)
    thermals.forward_euler(p_max)
    T_max = thermals.T_tank[-1]
    p_tank = n2o.p(T_max)

    return p_tank - p_measured

res = opt.root_scalar(func, bracket = bounds , method='bisect')
print(res)



