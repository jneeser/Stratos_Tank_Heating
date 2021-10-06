import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

class N2O:
    def __init__(self):
        self.T_c = 309.57   # K
        self.p_c = 7251.    # kPa
        self.rho_c = 452.   # kg/m3
        self.m = 174.       # kg
        self.V = 0.259      # m3
        self.hvap = 375e3   # J/kg      https://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=4

    def vap_mass_frac(self, T):
        # from SIII diptube calc 
        # linear extrapolation 
        # use this to check results

        temp = np.array([287,288,289,290,291,292,293,294,295,296,297,298,299,300,301])
        vap_frac_arr = np.array([0.0343,0.0356,0.037,0.0384,0.0399,0.0415,0.0431,0.0449,0.0468,0.0488,0.051,0.0533,0.0558,0.0585,0.0615])
        f = interpolate.interp1d(temp, vap_frac_arr, fill_value='extrapolate')

        vap_frac = f(T)

        return vap_frac


    def p(self, T):
        T_r = T / self.T_c
        a = 1 - T_r
        b = [-6.71893, 1.35966, -1.3779, -4.051]
        n = [1., 3./2., 5./2., 5.]
        p = 0
        for i in range(len(b)):
            p += b[i]*a**n[i]
        return self.p_c*np.exp(1/T_r*p)


    def rho_L(self, T):
        T_r = T / self.T_c
        a = 1 - T_r
        b = [1.72328, -0.83950, 0.51060, -0.10412]
        n = [1./3., 2./3., 3./3., 4./3.]
        rho = 0
        for i in range(len(b)):
            rho += b[i]*a**n[i]
        return self.rho_c*np.exp(rho)


    def rho_V(self, T):
        T_r = T / self.T_c
        a = 1/T_r-1
        b = [-1.00900, -6.28792, 7.50332, -7.90463, 0.629427]
        n = [1./3., 2./3., 3./3., 4./3., 5/3]
        rho = 0
        for i in range(len(b)):
            rho += b[i]*a**n[i]
        return self.rho_c*np.exp(rho)


    def h_L(self, T):
        T_r = T / self.T_c
        a = 1 - T_r
        b = [-200., 116.043, -917.225, 794.779, -589.587]
        n = [0, 1./3., 2./3., 3./3., 4./3.]
        h = 0
        for i in range(len(b)):
            h += b[i]*a**n[i]
        return h


    def h_V(self, T):
        T_r = T / self.T_c
        a = 1 - T_r
        b = [-200., 440.055, -459.701, 434.081, -485.338]
        n = [0, 1./3., 2./3., 3./3., 4./3.]
        h = 0
        for i in range(len(b)):
            h += b[i]*a**n[i]
        return h


    def cp_L(self, T):
        T_r = T / self.T_c
        a = 1 - T_r
        b = [2.49973, 0.023454, -3.80136, 13.0945, -14.5180]
        n = [0, -1., 1, 2, 3]
        cp = 1
        for i in range(1,len(b)):
            cp += b[i]*a**n[i]
        return b[0]*cp


    def cp_V(self, T):
        T_r = T / self.T_c
        a = 1 - T_r
        b = [132.632, 0.052187, -0.364923, -1.20233, 0.536141]
        n = [0, -2./3., -1./3., 1./3., 2./3.]
        cp = 1
        for i in range(1,len(b)):
            cp += b[i]*a**n[i]
        return b[0]*cp

    def alpha(self, T):
        return (self.m / self.V - self.rho_V(T)) / (self.rho_L(T) - self.rho_V(T))

    def m_V(self, T):
        return (1 - self.alpha(T)) * self.rho_V(T) * self.V

    def m_L(self, T):
        return self.alpha(T) * self.rho_L(T) * self.V

    def Q_vap(self, T, T_last): # Not the fastest but quite compact
        return (self.h_V(T)-self.h_L(T) - (self.p(T_last) - self.p(T))/self.rho_V(T_last)) * (self.m_L(T_last)-self.m_L(T))


if __name__ == '__main__':
    n2o = N2O()
    T = np.linspace(288, 309, 100)
    vap_frac = n2o.vap_mass_frac(T)

    plt.plot(T, vap_frac)
    plt.xlabel('temp [K]')
    plt.ylabel('vapor fraction [kg/kg]')
    plt.show()