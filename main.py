import numpy as np
import time

#import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#from astropy.time import Time
#from astroquery.jplhorizons import Horizons      

class TwoBodyEngine:
    
    def __init__(self, mu, step, timeEnd):

        # General Variables
        self.mu = mu                    # Gravitational Parameter for the central object ()
        self.step = step                # Time step
        self.timeStart = 0.0            # Starting time
        self.timeEnd = timeEnd          # Final time
        self.numberOfLoops = None       # Number of loops
        #self.Rworking = None            # Current Location Vector
        #self.Vworking = None            # Current Velocity Vector

        # Flags
        self.solarRadiationPressureFlag = 0         # Solar radiation pressure (on(1), off(0))
        self.thirdPlanetPerturbationFlag = 0        # Third planet perturbation (on(1), off(0))
        

        # Solar Radiation Pressure (Curtis page 695)
        self.nu = 1.0                   # Object Always in the sun
        self.solarCosntant = 63.15e12   # Solar Constant W/Km2
        self.speedOfLight = 2.988e5     # Speed of light Km/s
        self.cr = 1.3                   # 1 for black, 2 for 100% reflective
        self.absorbingArea = 1.0e-6     # Abdsorbing area, Pi*R**2 for connonbell modal

    def twoBodyEquations(self, R, V):

        # Two Body Equations
        RDot = np.array(V)
        RDotDot = np.array([-(self.mu / pow(np.linalg.norm(R), 3)) * i for i in R])
        return RDot, RDotDot

    def twoBodyEquationsPer(self, R, V):

        # Two Body Equations With Perturbations
        RDot = np.array(V)
        Perturbations = self.solarRadiationPressureFlag * self.solarRadiationPressure(R) + self.thirdPlanetPerturbationFlag * self.thirdPlanetPerturbation(R, )#TODO
        RDotDot = np.array([-(self.mu / pow(np.linalg.norm(R), 3.0)) * i for i in R]) + Perturbations
        return RDot, RDotDot

    def solarRadiationPressure(self, RHeliocentric):

        # Solar Radiation Pressure
        # DOES IT NEED -1 behind nu?
        accelerationSRP = self.nu * (self.solarCosntant * pow(696.0 / np.linalg.norm(RHeliocentric), 2.0)) * self.cr * (1.0 / self.speedOfLight) * self.absorbingArea * (RHeliocentric / np.linalg.norm(RHeliocentric))
        return accelerationSRP

    def thirdPlanetPerturbation(self, RHeliocentric, RThirdPlanet):
        # Third Planet Perturbations:
        RRelative = RThirdPlanet - RHeliocentric
        accelerationTPP = MUp * ((RRelative) / pow(np.linalg.norm(RRelative), 3) - (RThirdPlanet) / pow(np.linalg.norm(RThirdPlanet), 3))#TODO
        return accelerationTPP

    def determineNumberOfSteps(self):

        # Determine Number of Steps
        self.numberOfLoops = int((self.timeEnd - self.timeStart) / self.step)

    def solverRK4Fixed(self, RInitial, VInitial):

        # 4th Order Rungee-Kutta With Fixed Steps
        self.determineNumberOfSteps()
        currentTime = self.timeStart
        for i in range(self.numberOfLoops):
            print(i+1,"/" , self.numberOfLoops)
            k1_Rd, k1_Rdd = self.step * self.twoBodyEquations(RInitial, VInitial)
            k2_Rd, k2_Rdd = self.step * self.twoBodyEquations(RInitial + (0.5 * k1_Rd), VInitial + (0.5 * k1_Rd))
            k3_Rd, k3_Rdd = self.step * self.twoBodyEquations(RInitial + 0.5 * k2_Rd, VInitial + 0.5 * k2_Rdd)
            k4_Rd, k4_Rdd = self.step * self.twoBodyEquations(RInitial + k3_Rd, VInitial + k3_Rdd)
            Rnew = (1.0 / 6.0) * (k1_Rd + 2.0 * k2_Rd + 2.0 * k3_Rd + k4_Rd) + RInitial
            Vnew = (1.0 / 6.0) * (k1_Rdd + 2.0 * k2_Rdd + 2.0 * k3_Rdd + k4_Rdd) + VInitial
            RInitial = Rnew
            VInitial = Vnew
            currentTime = currentTime + self.step
            print("R:", Rnew)
            print("V:", Vnew)
            print("T:", currentTime, "Sec", "\n")

class Planet:

    def __init__(self, name, mu):

        # Object Planet
        self.name = name                # Name of the planet
        self.mu = mu                    # Gravitational Parameter Km3/m2
        self.r = None                   # Current Location Vector
        self.v = None                   # Current Velocity Vector  

""" class SolarSystem:

    def __init__(self, name, mu):

        # Solar System """
         

t0 = time.time()
print("\n")

# Defining Solar System
SS = TwoBodyEngine(0.5, 1, 3)
SS.solverRK4Fixed([1,2,3],[1,2,3])

tf = (time.time() - t0)
print("Elapsed Time:", tf, "Sec")
