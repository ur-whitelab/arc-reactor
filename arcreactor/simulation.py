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
        c = [0,0,0]
        '''Calculates all the kinetics of the reactors'''
        c0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        tau = self.reactor_volumes[0] / self.volumetric_feed_rates[0]
        k = 0.01
        t_int = np.linspace(0, 3600, 3600000)  # simulate for one hour

        def rxn_d(c, t, k=0.01):
            return -k * c

        ca1 = odeint(rxn_d, c0[0], t_int)
        x1 = 1 - ca1 / c0[0]    #Conversion
        x = int(simulation_state.time * 1000 - self.start_time * 1000)
        # The following steps for calculating reactor concentrations can be looped later
        c[0] = np.asarray([ca1[x], ca1[x], c0[0] - ca1[x], c0[0] - ca1[x]])
        mf1 = c[0] / np.sum(c[0])
        ca2 = ca1 / (1 + tau * k)
        c[1] = np.asarray([ca2[x], ca2[x], c0[0] - ca2[x], c0[0] - ca2[x]])
        mf2 = c[1] / np.sum(c[1])
        ca3 = ca2 / (1 + tau * k)
        c[2] = np.asarray([ca3[x], ca3[x], c0[0] - ca3[x], c0[0] - ca3[x]])
        mf3 = c[2] / np.sum(c[2])
        temp = 25
        pressure =  1                      #Specific to this reaction

        simulation_state.time = int(time.time())
        for i in range(3):
            rxr = ReactionKinetics()
            rxr.temperature = temp
            rxr.pressure = pressure
            rxr.mole_fraction = c[i]
            simulation_state.kinetics.add(rxr) 

        return simulation_state



