import numpy as np
import asyncio
import datetime as dt
import time
from scipy.integrate import odeint
from .protobufs.reactors_pb2 import ReactorSystem
from .protobufs.kinetics_pb2 import SystemKinetics

# Hydrolysis of ester - A psuedo first order chemical reaction
# All reactors are equally sized and operate at the same temperature
# C2H5COOCH3 + H2O --> CH3COOH + C2H5OH

class Simulation:

    def __init__(self):
        #For demo purposes, the values are fixed
        self.reactor_number = 3
        self.reactor_volumes = np.asarray([2,2,2])
        self.volumetric_feed_rates = np.asarray([0.5,0.5])
        self.molar_feed_rate = np.asarray([1,1])      

    def calculate(self):
        '''Calculates all the kinetics of the reactors'''
        c0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        tau = reactor_volumes[0] / feed_volumetric_rates[0]
        t_int = np.linspace(0, 3600, 3600000)  # simulate for one hour

        def rxn_d(c, t, k=0.01):
            return -k * c
        
        ca1 = odeint(rxn_d, c0[0], t_int)
        x1 = 1 - ca1 / c0[0]    #Conversion
        x = int(time.time() * 1000 - self.start_time * 1000)

        # The following steps for calculating reactor concentrations can be looped later
        c1 = np.asarray([ca1[x], ca1[x], c0[0] - ca1[x], c0[0] - ca1[x]])
        mf1 = c1 / np.sum(c1)
        ca2 = ca1 / (1 + tau * k)
        c2 = np.asarray([ca2[x], ca2[x], c0[0] - ca2[x], c0[0] - ca2[x]])
        mf2 = c2 / np.sum(c2)
        ca3 = ca2 / (1 + tau * k)
        c3 = np.asarray([ca3[x], ca3[x], c0[0] - ca3[x], c0[0] - ca3[x]])
        mf3 = c3 / np.sum(c3)
        temp = 'room temperature'
        pressure =  '1atm'                      #Specific to this reaction
        return c1, c2, c3, temp, pressure       #Concentration with time


