import numpy as np
import matplotlib.pyplot as plt
from iapws import IAPWS97
from scipy.optimize import newton

### All distances are in m
### All units of energy are in W or J
### All temperatures are in C
### Pressure is expected in Pa

class Reactor:

    def __init__(self, name, Hi, Dr, Pitch, 
                 Df, Gap, kg, kf, kc, G, q0, 
                 P0, Tf0, dz, mode='PWR'):

        self.name = name
        self.Hi = Hi
        self.Dr = Dr
        self.Pitch = Pitch
        self.Df = Df
        self.Gap = Gap
        self.kg = kg
        self.kf = kf
        self.kc = kc
        self.G = G
        self.q0 = q0
        self.P0 = P0 # Pa
        self.Tf0 = Tf0
        self.dz = dz
        self.mode = mode

        self.Pw = np.pi * self.Dr
        self.Ph = np.pi * self.Dr
        self.Af = self.Pitch**2 - (0.25 * np.pi * (self.Dr**2))
        self.Dh = 4 * self.Af / self.Ph #hydraulic diameter

        self.R1 = self.Dr / 2 # Clad outer boundary
        self.R3 = self.Df / 2 # Gap : Fuel boundary
        self.R2 = self.R3 + self.Gap # Clad : Gap boundary

        self.State_Initial = IAPWS97(P=self.P0/1e6, T=self.Tf0+273.15) # MPa, K
        self.Tsat = IAPWS97(P=self.P0/1e6, x=0).T - 273.15   # C
        self.cp = self.State_Initial.cp * 1e3 # J/kg·K
        self.kfluid = self.State_Initial.k # W/m·K
        self.Pr = self.mu * self.cp / self.kfluid
        self.Re = self.G * self.Dh / self.mu

        self.z = np.arange(0, self.Hi + (self.dz / 2), self.dz)

    def q1(self, z): # Linear Heat Generation in fuel
        return self.q0 * np.sin(np.pi * z / self.Hi)

    def q3(self, z): # Volumetric Heat Generation in fuel
        A_fuel = 0.25 * self.Df**2 * np.pi
        ret = self.q1(z) / A_fuel
        return ret

    def q2(self, z):
        return self.q1(z) / self.Ph

    def mu(self,z):
        return ((self.X(z)/IAPWS97(P=self.Pressure(z)/1e6,x=1).mu) + ((1-self.X(z)) / IAPWS97(P=self.Pressure(z)/1e6,x=0).mu))**(-1) # TODO check units

    #########################################
    #  Pressure Drop
    #########################################

    def Pressure(self,z):
        if self.mode == 'PWR':
            f = self.f1P()
            acc = 0
        if self.mode == 'BWR' :
            acc = 0 #TODO
        DP = f * ((self.G**2 * self.Pw) / (8 * self.rho * self.Af)) * z
        DP = DP + (self.rho * 9.81 * z) + acc
        return self.P0 - DP

    def f1P(self):
        return ((0.790 * np.log(self.Re)) - 1.64)**(-2)

    def f2P(self):
        return 0 # TODO

    ###########################
    #  Quality and Properties
    ###########################

    def X(self,z):
        if self.Xe(z) < 0:
            return 0
        if (self.Xe(z) >= 0) or (self.Xe(z) <= 1):
            return self.Xe(z)
        if self.Xe(z) > 1:
            return 1

    def Xe(self,z):
        h0 = IAPWS97(P=self.P0/1e6, T=self.Tf0 + 273.15).h
        sat = IAPWS97(P=self.P0/1e6, x=0)
        hfg = sat.Hvap
        Xe_in = (h0 - sat.h) / hfg

        A = (self.q0 / (self.Af * self.G * hfg)) * (self.Hi / np.pi)
        X = A * (1- np.cos(np.pi * z / self.Hi)) + Xe_in
        return X
             
    def rho_m(self,z):
        r = 1/(IAPWS97(P=self.Pressure(z),x=0).v + IAPWS97(P=self.Pressure(z),x=1).v* self.X(z))
        return r

    def T2P(self,z):
        T_wall = np.zeros_like(self.z)
        T_bulk = np.zeros_like(self.z)
        T_sat = np.zeros_like(self.z)
        Xe_arr = np.zeros_like(self.z)
        X_arr = np.zeros_like(self.z)
        P_arr = np.zeros_like(self.z)

        for i in len(T_wall):
            if i == 0: 
                T_bulk[i] = self.Tf0
                Xe_arr[i] = self.Xe(0)
                X_arr[i] = 0
                P_arr[i] = self.P0

            state = IAPWS97(P=self.Pressure(z)/1e6, T=self.T_bulk[i])
            sat = IAPWS97(P=self.Pressure(z)/1e6, x=0)

            if Xe_arr[i] < 0:
                hfg = IAPWS97(P=P_arr[i],x=1).h - IAPWS97(P=P_arr[i],x=0).h
                Tsat = IAPWS97(P=P_arr[i],x=0).T - 273.15
                TF = (Xe_arr[i] * hfg / self.cp) + Tsat
                
                T_bulk[i] = TF
                T_sat[i] = Tsat

            if TF[i-1] == T_sat[i-1]:
                return Tsat #TSAT


    ###########################
    #  Heat Convection in fluid
    ###########################

    def T1P(self,z):

        if self.mode == 'PWR' :
            C7 = self.Tf0
            C8 = (self.q0 * self.Hi) / (self.G * self.Af * self.cp * np.pi)
            T = C7 + C8 * (1 - np.cos(np.pi * z / self.Hi))

        return T

    def h(self):
        if self.mode == 'PWR':
            Nu = 0.023 * (self.Re**0.8) * (self.Pr**0.4)
            h = (self.kfluid * Nu) / self.Dh
        if self.mode == 'BWR':
            h = 0 #TODO
        return h
        
    #########################################
    #  Heat Generation and Conduction In Rods
    #########################################

    def T_Rod(self,r,z):
        q3 = self.q3(z)

        if self.mode == 'PWR':
            T_inf = self.T1P(z)
        if self.mode == 'BWR':
            T_inf = self.T2P(z)

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


Option_1 = Reactor(
    name='PWR Option 1',
    Hi=4.1, # m
    Dr=0.94e-2, # m
    Pitch=1.25e-2, # m
    Df=0.81e-2, # m
    Gap=0.005e-2, # m
    kg=0.25, # W/mC
    kf=3.6, # W/mC
    kc=21.5, # W/mC
    G= 3800, # kg/m2 s
    q0=380 * 1e2, # W/m
    P0=15e6, # Pa
    Tf0=277, # C
    dz = 0.1, # m 
    mode= 'PWR' )