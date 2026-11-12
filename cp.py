import numpy as np
import matplotlib.pyplot as plt
from iapws import IAPWS97
from scipy.optimize import newton

### All distances are in m
### All units of energy are in W or J
### All temperatures are in C
### Pressure is expected in Pa

class BWR:

    def __init__(self, name, Hi, Dr, Pitch, 
                 Df, Gap, kg, kf, kc, G, q0, 
                 P0, Tf0, dz):

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

        self.Pw = np.pi * self.Dr
        self.Ph = np.pi * self.Dr
        self.Af = self.Pitch**2 - (0.25 * np.pi * (self.Dr**2))
        self.Dh = 4 * self.Af / self.Ph #hydraulic diameter

        self.R1 = self.Dr / 2 # Clad outer boundary
        self.R3 = self.Df / 2 # Gap : Fuel boundary
        self.R2 = self.R3 + self.Gap # Clad : Gap boundary

        self.State_Initial = IAPWS97(P=self.P0/1e6, T=self.Tf0+273.15) # MPa, K
#        self.Tsat = IAPWS97(P=self.P0/1e6, x=0).T - 273.15   # C
        self.cp = self.State_Initial.cp * 1e3 # J/kg K
        self.kfluid = self.State_Initial.k # W/m K
        self.Pr = self.mu * self.cp / self.kfluid
        self.Re = self.G * self.Dh / self.mu
        self.mu = self.State_Initial.mu
        self.z = np.arange(0, self.Hi + (self.dz / 2), self.dz)
        self.hf_sat = self.h_sat()[0]
        self.hg_sat = self.h_sat()[1]
        self.hfg_sat = self.h_sat()[2]

    def q1(self, z): # Linear Heat Generation in fuel
        return self.q0 * np.sin(np.pi * z / self.Hi)

    def q3(self, z): # Volumetric Heat Generation in fuel
        A_fuel = 0.25 * self.Df**2 * np.pi
        ret = self.q1(z) / A_fuel
        return ret

    def q2(self, z): # Heat flux
        return self.q1(z) / self.Ph

    def h_sat(self): # Saturation entapies
        hf = IAPWS97(P=self.P0/1e6,x=0).h * 1e3
        hg = IAPWS97(P=self.P0/1e6,x=1).h * 1e3
        hfg = (hg - hf) * 1e3 # J/kg K
        return [hf, hg, hfg]


    #########################################
    #  Pressure Drop
    #########################################

    def Pressure(self,z):
        return 0

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

    def Xe(self, z):
        Xe_in = (self.cp * self.Tf0) * 1e3 / self.hfg_sat

        A = (self.q0 / (self.Af * self.G * self.hfg_sat)) * (self.Hi / np.pi)
        return A * (1 - np.cos(np.pi * z / self.Hi)) + Xe_in
             
    # def rho_m(self,z):
    #     vf = IAPWS97(P=self.Pressure(z),x=0).v
    #     vg = IAPWS97(P=self.Pressure(z),x=1).v
    #     Xe = self.Xe(z)

    #     if Xe <= 0:
    #         X = 0
    #     if 0 < Xe and Xe <= 1:
    #         X = Xe
    #     if Xe > 1:
    #         X = 1

    #     r = (X * vg + (1-X)*vf)**(-1)
    #     return r

    def T2P(self,z):
        T_wall = np.zeros_like(self.z)
        T_bulk = np.zeros_like(self.z)
        T_sat = np.zeros_like(self.z)
        Xe_arr = np.zeros_like(self.z)
        X_arr = np.zeros_like(self.z)
        P_arr = np.zeros_like(self.z)

        for i in len(T_wall):
            if i == 0:  # If you are on the first step, initialize
                T_bulk[i] = self.Tf0
                Xe_arr[i] = self.Xe(0)
                X_arr[i] = 0
                P_arr[i] = self.P0

            state = IAPWS97(P=self.Pressure(z)/1e6, T=self.T_bulk[i])
            sat = IAPWS97(P=self.Pressure(z)/1e6, x=0)

            if Xe_arr[i-1] < 0: # If your quality is bellow 0, you are subcooled
                hfg = IAPWS97(P=P_arr[i-1],x=1).h - IAPWS97(P=P_arr[i-1],x=0).h
                Tsat = IAPWS97(P=P_arr[i-1],x=0).T - 273.15
                TF = (Xe_arr[i-1] * hfg / self.cp) + Tsat
                
                T_bulk[i] = TF
                T_sat[i] = Tsat
                Xe_arr[i] = self.Xe(z=self.z[i])

                #TODO pressure drop part for P_arr
                P_arr[i] = 0 #TODO

            if (Xe_arr[i-1] >= 0) and (Xe_arr[i-1] >= 0):
                return 0 #TODO

    def Pressure_solve(self,z,P):
        Xe = self.Xe(z)
        # P will be like P_arr[i-1]

        if Xe < 0: # single phase fluid
            rho_f = IAPWS97(P=P/1e6,x=0).rho
            f = ((0.790 * np.log(self.Re)) - 1.64)**(-2)

            frictional = f * ((self.G**2 * self.Pw) / (8 * rho_f * self.Af)) * z
            gravitational = (rho_f * 9.81 * z) 
            DP = frictional + gravitational
            return self.P0 - DP

        if (Xe >= 0) and (Xe >= 1): # double phase 
            v_f = IAPWS97(P=P/1e6,x=0).v
            v_g = IAPWS97(P=P/1e6,x=1).v

            rho_m = ((Xe * v_g) + (1-Xe)*v_f)**(-1)



        if Xe > 1: # Superheated vapor 
            rho_f = IAPWS97(P=P/1e6,x=1).rho
            f = 0 #TODO find a correlation for fluid 

            frictional = f * ((self.G**2 * self.Pw) / (8 * rho_f * self.Af)) * z
            gravitational = (rho_f * 9.81 * z) 
            DP = frictional + gravitational
            return self.P0 - DP



    ###########################
    #  Heat Convection in fluid
    ###########################

    def Liu(self,z,P):
        Xe = self.Xe(z)
        state = IAPWS97(P=P/1e6,x=Xe)
        rho_f = IAPWS97(P=P/1e6,x=0).rho
        rho_g = IAPWS97(P=P/1e6,x=1).rho
        Pr = state.mu * state.cp / (state.k)

        F = (1 + (self.Xe(z) * Pr * ((rho_f / rho_g) -1)))**(0.35)
        Re = self.G * self.Dh / state.mu
        S = (1 + (0.055*F**0.1) * (Re**0.16))**(-1)

        Nu = 0.023 * (Re**0.8) * (Pr**0.4)
        h1p = (state.k * Nu) / self.Dh

        ## h2p correlation

    # def LiuCoeff(self,z,P):
    #     Xe = self.Xe(z)
    #     state = IAPWS97(P=P/1e6,x=Xe)
    #     rho_f = IAPWS97(P=P/1e6,x=0).rho
    #     rho_g = IAPWS97(P=P/1e6,x=1).rho
    #     Pr = state.mu * state.cp / (state.k)

    #     F = (1 + (self.Xe(z) * Pr * ((rho_f / rho_g) -1)))**(0.35)
    #     Re = self.G * self.Dh / state.mu
    #     S = (1 + (0.055*F**0.1) * (Re**0.16))**(-1)

    #     return (F,S)
        
    #########################################
    #  Heat Generation and Conduction In Rods
    #########################################

    def T_Rod(self,r,z,T_bulk,h):
        # Now accepts Tbulk from Tbulk_arr[i]
        # Accepts new h from liu-winterton 
        q3 = self.q3(z) 

        C1 = -(q3 * self.R3**2) / (2 * self.kc)
        C2 = ((q3 * self.R3**2) / (h * 2 * self.kc))
        C2 = C2 * ((self.kc / self.R1) + h * np.log(self.R1)) + T_bulk
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
