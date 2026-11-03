import numpy as np
import matplotlib.pyplot as plt

### All distances are in m
### All units of energy are in W or J


# saturation temp of water at 15MPA 
# https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab3.pdf

class PWR:

    def __init__(self,name,Hi,Dr,P,Df,G,kg,kf,kc,q0,P0,Tf0):

        self.name = 'PWR'
        self.Hi = 4.1 #m
        self.Dr = 0.94e-2 #m
        self.P = 1.25e-2 #m
        self.Df = 0.81e-2 #m
        self.G = 0.005e-2 #m
        self.kg = 0.25 #W/mC
        self.kf = 3.6 #W/mC
        self.kc = 21.5 #W/mC
        self.q0 = 380 * 1e2 #W/m
        self.P0 = 15e6 #Pa 
        self.Tf0 = 277 #C

        self.Tsat = 342.155 #C

    def q1(self, z):
        ret = self.q0 * np.sin(np.pi * z / self.Hi)
        return ret

    

    