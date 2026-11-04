import numpy as np
import matplotlib.pyplot as plt

### All distances are in m
### All units of energy are in W or J
### All temperatures are in C

# Saturation temp of water at 15MPA 
# https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab3.pdf


PWR_dimensions = ('Hi' = 4.1 , 
                  'Dr' = 0.94e-2 )

class PWR:

    def __init__(self, name, Hi, Dr, P, Df, G, kg, kf, kc, q0, P0, Tf0):

        self.name = name
        self.Hi = Hi
        self.Dr = Dr
        self.P = P
        self.Df = Df
        self.G = G
        self.kg = kg
        self.kf = kf
        self.kc = kc
        self.q0 = q0
        self.P0 = P0
        self.Tf0 = Tf0

        self.Tsat = 342.155

    def q1(self, z):
        return self.q0 * np.sin(np.pi * z / self.Hi)

    


PWR1 = PWR(
    name='PWR 1',
    Hi=4.1, # m
    Dr=0.94e-2, # m
    P=1.25e-2, # m
    Df=0.81e-2, # m
    G=0.005e-2, # m
    kg=0.25, # W/mC
    kf=3.6, # W/mC
    kc=21.5, # W/mC
    q0=380 * 1e2, # W/m
    P0=15e6, # Pa
    Tf0=277 ) # C 


