import numpy as np
import matplotlib.pyplot as plt

### All distances are in m
### All units of energy are in W or J
### All temperatures are in C

# Saturation temp of water at 15MPA 
# https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab3.pdf

class PWR:

    def __init__(self, name, Hi, Dr, Pitch, Df, G, kg, kf, kc, q0, P0, Tf0, dz):

        self.name = name
        self.Hi = Hi
        self.Dr = Dr
        self.Pitch = Pitch
        self.Df = Df
        self.G = G
        self.kg = kg
        self.kf = kf
        self.kc = kc
        self.q0 = q0
        self.P0 = P0
        self.Tf0 = Tf0
        self.dz = dz

        self.Tsat = 342.155

        self.Temp = np.arange(0,self.Hi,dz)
        self.Temp[0] = self.Tf0

        self.P = np.arrange(0,self.Hi,dz)
        self.P[0] = self.P0

        self.Xe = np.arange(0,self.Hi,dz)

        self.Pw = np.pi * self.Dr
        self.Ph = np.pi * self.Dr
        self.Af = self.Pitch**2 - (0.25 * np.pi * (self.Dr**2))
        self.Dh = 4 * self.Af / self.Ph

    def q1(self, z): # Linear Heat Generation in fuel
        return self.q0 * np.sin(np.pi * z / self.Hi)

    def q3(self, z): # Volumetric Heat Generation in fuel
        A_fuel = 0.25 * self.Df**2 * np.pi
        ret = self.q1(z) / A_fuel
        return ret


    ##########
    #  Heat Generation In Rods
    ##########
    
    # [3]

    def T_clad(self,z,r):
        R = self.Dr / 2
        C1 = 
        T = C1 * np.log(r) + C2


    # [2]


    # [1] Fuel

    # def T_fuel(self,z,r):
    #     T = -(self.q3(z) * r**2) / (4 * self.kf)
    #     C = 
    #     return T 
    
    # [2] Rod

#    def T_gap(self,z):




    # [3] Clad



    # T end




# PWR1 = PWR(
#     name='PWR 1',
#     Hi=4.1, # m
#     Dr=0.94e-2, # m
#     P=1.25e-2, # m
#     Df=0.81e-2, # m
#     G=0.005e-2, # m
#     kg=0.25, # W/mC
#     kf=3.6, # W/mC
#     kc=21.5, # W/mC
#     q0=380 * 1e2, # W/m
#     P0=15e6, # Pa
#     Tf0=277 ) # C 


