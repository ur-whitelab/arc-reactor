import numpy as np
import asyncio
import datetime as dt
import time
from scipy.integrate import odeint
from .protobufs.kinetics_pb2 import *


# Hydrolysis of ester - A psuedo first order chemical reaction
# All reactors are equally sized and operate at the same temperature
# C2H5COOCH3 + H2O --> CH3COOH + C2H5OH

class Simulation:

    def __init__(self, start_time):
        #For demo purposes, the values are fixed
        self.reactor_number = 3
        self.reactor_volumes = np.asarray([2,2,2])
        self.volumetric_feed_rates = np.asarray([0.5,0.5])
        self.molar_feed_rate = np.asarray([1,1])
        self.start_time = start_time

    def calculate(self, simulation_state):
        c = []
        '''Calculates all the kinetics of the reactors'''
        c0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        tau = self.reactor_volumes[0] / self.volumetric_feed_rates[0]
        k = 0.01
        t_int = np.linspace(0, 3600, 3600000)  # simulate for one hour -- compensate for frame-by-frame updates

        def rxn_d(c, t, k=0.01):
            return -k * c

        ca1 = odeint(rxn_d, c0[0], t_int)
        ca2 = ca1 / (1 + tau * k)
        ca3 = ca2 / (1 + tau * k)
        x1 = 1 - ca1 / c0[0]    #Conversion
        x = int(simulation_state.time * 1000 - self.start_time * 1000)
        # The following steps for calculating reactor concentrations can be looped later

        val1 = float(ca1[x])
        val2 = float(ca2[x])
        val3 = float(ca3[x])
        c.append([val1, val1, c0[0] - val1, c0[0] - val1])
        mf1 = c[0] / sum(c[0])
        c.append([val2, val2, c0[0] - val2, c0[0] - val2])
        mf2 = c[1] / sum(c[1])
        c.append([val3, val3, c0[0] - val3, c0[0] - val3])
        mf3 = c[2] / sum(c[2])
        temp = 25
        pressure =  1                      #Specific to this reaction

        simulation_state.time = int(time.time())
        for i in range(len(c)):
            rxr = ReactorKinetics()
            rxr.temperature = temp
            rxr.pressure = pressure
            for j in range(len(c[i])):
                rxr.mole_fraction.append(c[i][j])
            simulation_state.kinetics.extend([rxr])
            #simulation_state.kinetics[i] = rxr
            print(simulation_state.kinetics)

        return simulation_state



