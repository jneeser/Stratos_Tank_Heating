import scipy.optimize as opt

from run import k_carbon, k_alu, t_carbon, t_alu, A_tank, p_max, emissivity, lamp_temp 
from tank_temp import Tank, Lamp, Thermals
from n2o_class import N2O

T0 = 288                # K
time = 10*60            # s
p_measured = 5700       # kPa
dt = 10                 # s
SF = 0.01               # intial guess for shape factor
bounds = [0.00001, 0.01]

# find shape fator that corresponds with these measured values 

tank = Tank(k_carbon, t_carbon, k_alu, t_alu, A_tank)
lamp = Lamp(emissivity, lamp_temp, SF)
n2o = N2O()             

def func(shape_factor):
    lamp = Lamp(emissivity, lamp_temp, shape_factor)
    thermals = Thermals(tank, lamp, n2o, time, dt, T0, evap=False)
    thermals.forward_euler(p_max)
    T_max = thermals.T_tank[-1]
    p_tank = n2o.p(T_max)

    return p_tank - p_measured

res = opt.root_scalar(func, bracket = bounds , method='bisect')
print(res)


