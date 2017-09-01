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
        if(len(simulation_state.kinetics) == 0):
            return simulation_state
        conc = []
        '''calculates all the kinetics of the reactors'''
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        tau = self.reactor_volumes[0] / self.volumetric_feed_rates[0]
        k = 0.01
        t_int = np.linspace(0, 3600, 3600000)  # simulate for one hour -- compensate for frame-by-frame updates

        def rxn_d(conc, t, k=0.01):
            return -k * conc

        x = int(simulation_state.time * 1000 - self.start_time * 1000)
        conc_limiting, vals = [], []
        for i in range(len(simulation_state.kinetics)):
            if i == 0:
                conc_limiting.append(odeint(rxn_d, conc0[0], t_int))
            else:
                conc_limiting.append(conc_limiting[i-1] / (1 + tau * k))
            vals.append(float(conc_limiting[i][x]))
            conc.append([vals[i], vals[i], conc0[0] - vals[i], conc0[0] - vals[i]])
        # The following steps for calculating reactor concentrations can be looped later
        temp = 25
        pressure =  1                      #Specific to this reaction

        #simulation_state.time = int(time.time())
        for i in range(len(simulation_state.kinetics)):#conc is the list of lists of concentrations of reactor species. its length is the number of reactors.
            simulation_state.kinetics[i].temperature = temp
            simulation_state.kinetics[i].pressure = pressure
            while(len(simulation_state.kinetics[i].mole_fraction) < len(conc[i])):
                simulation_state.kinetics[i].mole_fraction.append(float(0))
            for j in range(len(conc[i])):
                simulation_state.kinetics[i].mole_fraction[j] = (conc[i][j])
        return simulation_state



