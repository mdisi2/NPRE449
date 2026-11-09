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

        ### Massive difference in BWR, will have to change due to pressure

        self.State_Initial = IAPWS97(P=self.P0/1e6, T=self.Tf0+273.15) # MPa, K
        self.Tsat = IAPWS97(P=self.P0/1e6, x=0).T - 273.15   # C
        self.mu = self.State_Initial.mu              # Pa *s
        self.cp = self.State_Initial.cp * 1e3        # J/kg·K
        self.rho = self.State_Initial.rho             # kg/m3
        self.kfluid = self.State_Initial.k               # W/m·K
        self.Pr = self.mu * self.cp / self.kfluid
        self.Re = self.G * self.Dh / self.mu

        ###

        self.z = np.arange(0, self.Hi, self.dz)
        self.P_arr = self.Pressure(self.z) # Pa
        self.Tsat_arr = np.array([IAPWS97(P=p/1e6, x=0).T - 273.15 for p in self.P_arr]) # C
        self.T_fluid_arr = np.zeros_like(self.z)
        self.T_fluid_arr[0] = self.Tf0

        self.Xe_arr = np.zeros_like(self.z)
        self.Xe_arr[0] = 0 # TODO

    def q1(self, z): # Linear Heat Generation in fuel
        return self.q0 * np.sin(np.pi * z / self.Hi)

    def q3(self, z): # Volumetric Heat Generation in fuel
        A_fuel = 0.25 * self.Df**2 * np.pi
        ret = self.q1(z) / A_fuel
        return ret

    def q2(self, z):
        return self.q1(z) / self.Ph

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

    def Pressure_Plot(self):

        plt.figure()
        plt.plot((self.P_arr)/1e6,self.z,label='Fluid Pressure',color='orange')
        plt.xlabel('Channel Pressure [MPa]')
        plt.ylabel('Channel Height [m]')
        plt.grid(alpha=0.5,color='gray',linestyle=":")
        plt.legend()
        plt.savefig(f'{self.name}_Pressure.png',dpi=600)

    ###########################
    #  Quality and Properties
    ###########################

    def T2P(self,z,i):
        TF = self.T_fluid_arr[i]
        Tsat = self.Tsat_arr[i]

        state = IAPWS97(P=self.Pressure(z)/1e6, T=self.T_fluid_arr[i])
        sat = IAPWS97(P=self.Pressure(z)/1e6, x=0)

        if TF < Tsat:
            X = 0
            Xe = (state.cp * TF) - (sat.cp * Tsat) / (sat.Hvap) 
            rho = ((1-X)* state.v + X*state.v)**(-1)
            T = (self.hfg(z) / state.cp) * Xe + Tsat 

            return T

        if TF == Tsat:
            X = Xe
            rho = ((1-X)* state.v + X*state.v)**(-1)
            return Tsat #TSAT
        
        if TF > Tsat:
            Xe = (state.cp * TF) - (sat.cp * Tsat) / (sat.Hvap)
            X = 1 
            rho = ((1-X)* state.v + X*state.v)**(-1)

            T = (self.hfg(z) / state.cp) * Xe + Tsat 

            return T

    def hf(self,z):
            state = IAPWS97(P=self.Pressure(z)/1e6,x=0) 
            return state.h * 1e3 # J/kg

    def hfg(self,z):
            state = IAPWS97(P=self.Pressure(z)/1e6, x=0) 
            return state.Hvap * 1e3 # K/kg

    def dPdz(self, z):
        dPdz_arr = np.gradient(self.P_arr, self.z)
        return np.interp(z, self.z, dPdz_arr)

    def dXedz(self, z):
        Xe = self.Xe(z)
        P = self.Pressure(z)/1e6 #MPa
        dP = 1e-3 # small step for delta p MPA
        dhf_dP = (IAPWS97(P=P+dP, x=0).h - IAPWS97(P=P-dP, x=0).h) * 1e3 / (2*dP*1e6)
        dhfg_dP = (IAPWS97(P=P+dP, x=1).Hvap - IAPWS97(P=P-dP, x=1).Hvap) * 1e3 / (2*dP*1e6)

        dP_dz = self.dPdz(z)
        dh_dz = self.q1(z) / (self.Af * self.G)
        return (dh_dz - (dhf_dP + Xe * dhfg_dP) * dP_dz) / self.hfg(z)

    def Xe0(self):
        h0 = IAPWS97(P=self.P0/1e6, T=self.Tf0 + 273.15).h
        sat = IAPWS97(P=self.P0/1e6, x=0)

        return (h0 - sat.h) / sat.Hvap

    #newton

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

    def Temperature_Plot(self):
        z = self.z

        plt.figure()

        # T = np.zeros_like(z)
        # for j, zj in enumerate(z):
        #     T[j] = self.T_Rod(r=0,z=zj)
        # plt.plot(T,z,label=f'Centerline',color='red')

        T = np.zeros_like(z)
        for j, zj in enumerate(z):
            T[j] = self.T_Rod(r=self.R1,z=zj)
        plt.plot(T,z,label=f'Clad Surface',color='orange')

        T = np.zeros_like(z)
        for j, zj in enumerate(z):
            T[j] = self.T_fluid(z=zj)
        plt.plot(T,z,label=f'Bulk Fluid',color='blue')

        plt.plot(self.Tsat_arr,z,label=f'Tsat',color='green',linestyle=':')
        plt.xlabel('Temperature [C]')
        plt.ylabel('Channel Height [m]')
        plt.grid(alpha=0.5,color='gray',linestyle=":")
        plt.legend()
        plt.savefig(f'{self.name}_Temp.png',dpi=600)

    def plotter(self): # returns all plots and questions
        self.Temperature_Plot()
        self.Pressure_Plot()
        print(f'{self.name} Plotted!')


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

Option_1.plotter()