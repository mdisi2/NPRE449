import numpy as np
import matplotlib.pyplot as plt
from iawps import IAPWS97

### All distances are in m
### All units of energy are in W or J
### All temperatures are in C

# Saturation temp of water at 15MPA 
# https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab3.pdf

class PWR:

    def __init__(self, name, Hi, Dr, Pitch, Df, Gap, kg, kf, kc, G, q0, P0, Tf0, dz):

        self.name = name
        self.Hi = Hi
        self.Dr = Dr
        self.Pitch = Pitch
        self.Df = Df
        self.Gap = Gap
        self.kg = kg
        self.kf = kf
        self.kc = kc
        self.q0 = q0
        self.P0 = P0
        self.Tf0 = Tf0
        self.dz = dz

        self.Tsat = 342.155 # C

        self.Temp_f = np.arange(0,self.Hi,dz)
        self.Temp_f[0] = self.Tf0

        self.P = np.arrange(0,self.Hi,dz)
        self.P[0] = self.P0

        self.Xe = np.arange(0,self.Hi,dz)

        self.Pw = np.pi * self.Dr
        self.Ph = np.pi * self.Dr
        self.Af = self.Pitch**2 - (0.25 * np.pi * (self.Dr**2))
        self.Dh = 4 * self.Af / self.Ph #hydraulic diameter


        self.R1 = self.Dr / 2 # Clad outer boundary
        self.R3 = self.Df / 2 # Gap : Fuel boundary
        self.R2 = self.R3 + self.Gap # Clad : Gap boundary

        ### Massive difference in BWR, will have to change due to pressure

        self.State_Initial = IAPWS97(P=self.P0/1e6, T=Tf0+273.15) #expects MPa and K
        self.mu = self.State_Initial.mu #constant viscosity Pa *s
        self.Re = self.Gap * self.Dh / self.mu
        self.cp = self.State_Initial.cp #constant cp kJ/kg·K
        self.rho = self.State_Initial.rho # constant density in single phase
        self.kfluid = self.State_Initial.ThCond() # constant thermal conductivity of water
        self.Pr = self.mu * self.cp / self.kfluid

    def q1(self, z): # Linear Heat Generation in fuel
        return self.q0 * np.sin(np.pi * z / self.Hi)

    def q3(self, z): # Volumetric Heat Generation in fuel
        A_fuel = 0.25 * self.Df**2 * np.pi
        ret = self.q1(z) / A_fuel
        return ret

    ###########################
    #  Heat Convection in fluid
    ###########################

    def T_fluid(self,z):
        C7 = self.Tf0
        C8 = (self.q0 * self.Hi) / (self.G * self.Af * self.cp * np.pi)
        T = C7 + C8 * (1 - np.cos(np.pi * z / self.Hi))
        return T

    def h(self,z):
        Nu = 0.023 * (self.Re**0.8) * (self.Pr**0.4)
        h = (self.kfluid * Nu) / self.Dh
        return h
        

    ###########################
    #  Heat Generation and Conduction In Rods
    ###########################

    def T_Rod(self,r,z):
        q3 = self.q3(z)
        T_inf = self.T_fluid(z)
        h = self.h(z)


        C1 = -(q3 * self.R3**2) / (2 * self.kc)
        C2 = ((q3 * self.R3**2) / (h * 2 * self.kc))
        C2 = C2 * ((self.kc / self.R1) + h * np.log(self.R1)) + T_inf
        C3 = -(q3 * self.R3**2) / (2 * self.kg)
        C4 = np.log(self.R2) * (C1 - C3) + C2

        C6 = C3*np.log(self.R3) + C4 + ((q3 * self.R3**2) / (4 * self.kf))

        if r <= self.R3: # in fuel
            T = ( (-q3 * r**2) / (4* self.kf) ) + C6

        if (r > self.R3) and (r <= self.R2): # in gap
            T = C3 * np.log(r) + C4

        if (r > self.R2) and (r <= self.R1): #in clad
            T = C1 * np.log(r) + C2

        return T


######
#   Conservation Equations
######

# def p_i(self,z) #solves for pressure at i of z
# def h_i(self,z) #solves for enthaply at i of z

# def Single_Phase(self,z):
#   nested equations T_clad(R=Diameter) \ T_gap(r) \ T_fuel(r)
#   Now can solve q''_clad into water coolant
#   single phase equations as a function of z
#   returns complete T_fluid and T_pressure array
#  two phase equations will return all that and T array 


# PWR1 = PWR(
#     name='PWR 1',
#     Hi=4.1, # m
#     Dr=0.94e-2, # m
#     P=1.25e-2, # m
#     Df=0.81e-2, # m
#     Gap=0.005e-2, # m
#     kg=0.25, # W/mC
#     kf=3.6, # W/mC
#     kc=21.5, # W/mC
#     G= 3800 # kg/m2 s
#     q0=380 * 1e2, # W/m
#     P0=15e6, # Pa
#     Tf0=277 ) # C