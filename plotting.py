import numpy as np
from matplotlib import pyplot as plt

def plot(pressure_list):
    data= np.genfromtxt('pressure.csv', delimiter=',', dtype=None, skip_header = 1)

    for i in range(1,len(pressure_list)):
        plt.plot(data[:,0]/60, data[:,i], label='P0 = '+str(pressure_list[i-1])+' kPa')

    plt.ylim([3400, 6100])
    plt.xlabel('time [min]')
    plt.ylabel('pressure [kPa]')
    plt.title('SIV Heating Time with 8 lamps')
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    
    p_range = np.arange(3500, 5500+200, 200)
    plot(p_range)