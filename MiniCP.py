import numpy as np
import matplotlib.pyplot as plt
from iapws import IAPWS97
from scipy.optimize import newton

### All distances are in m
### All units of energy are in W or J
### All temperatures are in C
### Pressure is expected in Pa

class PWR:

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
        self.Tsat = IAPWS97(P=self.P0/1e6, x=0).T - 273.15 # C
        self.mu = self.State_Initial.mu # Pa *s
        self.cp = self.State_Initial.cp * 1e3 # J/kg K
        self.rho = self.State_Initial.rho # kg/m3
        self.kfluid = self.State_Initial.k # W/m K
        self.Pr = self.mu * self.cp / self.kfluid
        self.Re = self.G * self.Dh / self.mu

        ###

        self.z = np.arange(0, self.Hi, self.dz)
        self.P_arr = self.Pressure(self.z) # Pa
        self.Tsat_arr = np.array([IAPWS97(P=p/1e6, x=0).T - 273.15 for p in self.P_arr]) # C

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
        f = self.f1P()
        DP = f * ((self.G**2 * self.Pw) / (8 * self.rho * self.Af)) * z
        DP = DP + (self.rho * 9.81 * z) 
        return self.P0 - DP

    def f1P(self):
        return ((0.790 * np.log(self.Re)) - 1.64)**(-2)

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

    def hfg(self,z):
        P = self.Pressure(z) /1e6
        hf = IAPWS97(P=P, x=0).h
        hg = IAPWS97(P=P, x=1).h
        return (hg - hf) * 1e3

    def Xe(self, z):
        h0 = IAPWS97(P=self.P0/1e6, T=self.Tf0 + 273.15).h
        sat = IAPWS97(P=self.P0/1e6, x=0)
        hfg0 = (IAPWS97(P=self.P0/1e6, x=1).h - sat.h) * 1e3   # J/kg
        Xe_in = (h0 - sat.h) * 1e3 / hfg0

        A = (self.q0 / (self.Af * self.G * hfg0)) * (self.Hi / np.pi)
        return A * (1 - np.cos(np.pi * z / self.Hi)) + Xe_in

    def Xe_Plotter(self):
        plt.figure()
        plt.plot([self.Xe(i) for i in self.z], self.z,
                label=r'Equilibrium Quality $\chi_E$', color='maroon')
        plt.xlabel(r'Equilibrium Quality $\chi_E$')
        plt.ylabel('Channel Height [m]')
        plt.grid(alpha=0.5,color='gray',linestyle=":")
        plt.legend()
        plt.savefig(f'{self.name}_Xe.png',dpi=600)

    
    ###########################
    #  Heat Convection in fluid
    ###########################

    def T_fluid(self,z):
        P = self.Pressure(z) / 1e6
        hfg = (IAPWS97(P=P,x=1).h - IAPWS97(P=P,x=0).h ) * 1e3
        Tsat = IAPWS97(P=P,x=0).T - 273.15
        TF = (self.Xe(z) * hfg / self.cp) + Tsat
        return TF

    def h(self,z):
        Nu = 0.023 * (self.Re**0.8) * (self.Pr**0.4)
        h = (self.kfluid * Nu) / self.Dh
        return h
        
    #########################################
    #  Heat Generation and Conduction In Rods
    #########################################

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

    def Temperature_Plot(self):
        z = self.z

        plt.figure()
        plt.xlabel('Temperature [C]')
        plt.ylabel('Channel Height [m]')
        plt.grid(alpha=0.5,color='gray',linestyle=":")

        T = np.zeros_like(z)
        for j, zj in enumerate(z):
            T[j] = self.T_Rod(r=self.R1,z=zj)
        plt.plot(T,z,label=f'Clad Surface',color='orange')

        T = np.zeros_like(z)
        for j, zj in enumerate(z):
            T[j] = self.T_fluid(z=zj)
        plt.plot(T,z,label=f'Bulk Fluid',color='blue')

        plt.plot(self.Tsat_arr,z,label=f'Tsat',color='green',linestyle=':')
       
        plt.legend()
        plt.savefig(f'{self.name}_Temp.png',dpi=600)

        T = np.zeros_like(z)
        for j, zj in enumerate(z):
            T[j] = self.T_Rod(r=0,z=zj)
        plt.plot(T,z,label=f'Centerline',color='red')
        
        plt.legend()
        plt.savefig(f'{self.name}_CenterlineTemp.png',dpi=600)


    def plotter(self): # returns all plots and questions
        self.Temperature_Plot()
        self.Pressure_Plot()
        self.Xe_Plotter()
        print(f'{self.name} Plotted!')
        plt.show()

Option_1 = PWR(
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
    dz = 0.1) # m

Option_1.plotter()

### Varying P

P = np.arange(15e6, 5e6, -0.1e6)

for i in P:
    reactor = PWR(
        name='PWR',
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
        P0=i, # Pa
        Tf0=277, # C
        dz = 0.1) # m

    xe = reactor.Xe(reactor.z)

    if np.max(xe) >= 0:
        break

plt.figure()
plt.plot(Option_1.Pressure(Option_1.z)/1e6, Option_1.z, label='Option 1 15MPA')
plt.plot(reactor.Pressure(reactor.z)/1e6, reactor.z, label=f'Critical {i/1e6:.3} MPa')
plt.xlabel('Pressure MPa')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Pressure Pressure.png',dpi=300)

plt.figure()
plt.plot([Option_1.T_fluid(zj) for zj in Option_1.z], Option_1.z, label='Option 1 bulk fluid')
plt.plot([reactor.T_fluid(zj) for zj in reactor.z], reactor.z, label=f'Critical {i/1e6:.3} MPa')
plt.axvline(x=Option_1.Tsat, color='gray', linestyle=':', label='Tsat Option 1')
plt.axvline(x=reactor.Tsat, color='k', linestyle=':', label='Tsat critical case')
plt.xlabel('Temp C')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Pressure bulk fluid.png',dpi=300)

plt.figure()
plt.plot([Option_1.T_Rod(0, zj) for zj in Option_1.z], Option_1.z, label='Option 1 centerline 15 MPa')
plt.plot([reactor.T_Rod(0, zj) for zj in reactor.z], reactor.z, label=f'Critical {i/1e6:.3} MPa Centerline')
plt.xlabel('Temp C')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Pressure centerline.png',dpi=300)

#pitch

Pitches = np.arange(1.25e-2, 0.95e-2, -0.005e-2)

for pitch in Pitches:
    reactor_p = PWR(
        name='PWR',
        Hi=4.1, # m
        Dr=0.94e-2, # m
        Pitch=pitch, # m
        Df=0.81e-2, # m
        Gap=0.005e-2, # m
        kg=0.25, # W/mC
        kf=3.6, # W/mC
        kc=21.5, # W/mC
        G= 3800, # kg/m2 s
        q0=380 * 1e2, # W/m
        P0=15e6, # Pa
        Tf0=277, # C
        dz = 0.1) # m

    xe = reactor_p.Xe(reactor_p.z)

    if np.max(xe) >= 0:
        break

pitch_crit = pitch

plt.figure()
plt.plot(Option_1.Pressure(Option_1.z)/1e6, Option_1.z, label='Option 1 pitch' )
plt.plot(reactor_p.Pressure(reactor_p.z)/1e6, reactor_p.z, label=f'Critical Geometry {pitch*1e2:.3}')
plt.xlabel('Pressure MPa')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Pitch Pressure.png',dpi=300)

plt.figure()
plt.plot([Option_1.T_fluid(zj) for zj in Option_1.z], Option_1.z, label='Option 1 bulk fluid')
plt.plot([reactor_p.T_fluid(zj) for zj in reactor_p.z], reactor_p.z, label=f'Critical Geometry {pitch*1e2:.3}')
plt.axvline(x=Option_1.Tsat, color='gray', linestyle=':', label='Tsat')
plt.xlabel('Temp C')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Pitch bulk fluid.png',dpi=300)

plt.figure()
plt.plot([Option_1.T_Rod(0, zj) for zj in Option_1.z], Option_1.z, label='Option 1 centerline')
plt.plot([reactor_p.T_Rod(0, zj) for zj in reactor_p.z], reactor_p.z, label=f'Critical Geometry {pitch*1e2:.3} Centerline')
plt.axvline(x=Option_1.Tsat, color='gray', linestyle=':', label='Tsat')
plt.xlabel('Temp C')
plt.ylabel('Channel Height')
plt.legend()
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.savefig('Critical Pitch centerline.png',dpi=300)
plt.show()

Tins = np.arange(277, 342, 0.5)   # C, stays below Tsat at 15 MPa (~342 C)

for Tin in Tins:
    reactor_T = PWR(
        name='PWR',
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
        Tf0=Tin, # C
        dz = 0.1) # m

    xe = reactor_T.Xe(reactor_T.z)

    if np.max(xe) >= 0:
        break

Tin_crit = Tin

plt.figure()
plt.plot(Option_1.Pressure(Option_1.z)/1e6, Option_1.z, label='Option 1 Tin 277 C')
plt.plot(reactor_T.Pressure(reactor_T.z)/1e6, reactor_T.z, label=f'Critical Tin {Tin_crit:.1f} C')
plt.xlabel('Pressure MPa')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Tin Pressure.png',dpi=300)

plt.figure()
plt.plot([Option_1.T_fluid(zj) for zj in Option_1.z], Option_1.z, label='Option 1 bulk fluid')
plt.plot([reactor_T.T_fluid(zj) for zj in reactor_T.z], reactor_T.z, label=f'Critical Tin {Tin_crit:.1f} C')
plt.axvline(x=Option_1.Tsat, color='gray', linestyle=':', label='Tsat')
plt.xlabel('Temp C')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Tin bulk fluid.png',dpi=300)

plt.figure()
plt.plot([Option_1.T_Rod(0, zj) for zj in Option_1.z], Option_1.z, label='Option 1 centerline')
plt.plot([reactor_T.T_Rod(0, zj) for zj in reactor_T.z], reactor_T.z, label=f'Critical Tin {Tin_crit:.1f} C Centerline')
plt.axvline(x=Option_1.Tsat, color='gray', linestyle=':', label='Tsat')
plt.xlabel('Temp C')
plt.ylabel('Channel Height')
plt.grid(alpha=0.5,color='gray',linestyle=":")
plt.legend()
plt.savefig('Critical Tin centerline.png',dpi=300)
plt.show()