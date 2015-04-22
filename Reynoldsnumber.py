# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:19:05 2015
Reynolds number based on temperature
@author: Roeland
"""
import numpy as np

#list of inputs: mach, (gamma, R), T, (Sutherland), pressure

class flowproperties:
    def __init__( self, temperaturearray, totalpressure):
        self.temperaturearray = temperaturearray
        self.totalpressure = totalpressure
        self.gamma = 1.4
        self.R = 287.05 
        self.M = 7.5 #mach
        self.T = temperaturearray*(1+(self.gamma-1)/2*self.M**2)**-1
        self.a = np.sqrt(self.gamma*self.R*self.T) #speed of sound
        self.V = self.M*self.a #velocity
    
        #base values for sutherland
        C = 120
        T0 = 291.15
        mu0 = 0.00001827
        
        self.P = self.totalpressure*(1+(self.gamma-1)/2*self.M**2)**(-self.gamma/(self.gamma-1))
        self.mu = mu0*(T0+C)/(self.T+C)*(self.T/T0)**1.5 #sutherland #kinematic viscosity
        self.rho = self.P/(self.R*self.T) #density
        self.v = self.mu/self.rho #dynamic viscosity
        self.chord = 01.00#meter
        self.Re = self.V*self.chord/self.v
        self.reynolds = self.Re
       
