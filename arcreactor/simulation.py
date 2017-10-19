import numpy as np
import asyncio
import datetime as dt
import time
from scipy.integrate import odeint
#from .protobufs.kinetics_pb2 import *


# Hydrolysis of ester - A psuedo first order chemical reaction
# All reactors are equally sized and operate at the same temperature
# C2H5COOCH3 + H2O --> CH3COOH + C2H5OH

class Simulation:
    '''Controls simulation of objects'''
    def __init__(self, start_time):
        #For demo purposes, the values are fixed
        self.reactor_number = 0
        self.reactor_volumes = np.asarray([2,2,2])
        self.volumetric_feed_rates = np.asarray([0.5,0.5])
        self.molar_feed_rate = np.asarray([1,1])
        self.start_time = start_time

    def calculate(self, simulation_state):
        '''Function that does the actual simulation for number of objects specified by the protobuf

        Parameters
        ----------
        simulation_state : protobuf object

        Returns
        -------
        protobuf object
                        Protobuf that specifies object properties

        '''
        if(len(simulation_state.kinetics) == 0):
            return simulation_state
        if(self.reactor_number != len(simulation_state.kinetics)):  #reset simulation when no of reactors change
            self.reactor_number = len(simulation_state.kinetics)
            self.start_time = simulation_state.time
        conc = []
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        tau = self.reactor_volumes[0] / self.volumetric_feed_rates[0]
        k = 0.01
        t_int = np.linspace(0, 3600, 3600*25)  # simulate for one hour -- compensation for frame-by-frame updates(~25fps)

        def rxn_d(conc, t, k=0.01):
            '''Calculating the rate of reaction'''
            return -k * conc

        x = simulation_state.time - self.start_time
        conc_limiting, vals = [], []
        for i in range(len(simulation_state.kinetics)):
            if i == 0:
                conc_limiting.append(odeint(rxn_d, conc0[0], t_int))
            else:
                conc_limiting.append(conc_limiting[i-1] / (1 + tau * k))
            vals.append(float(conc_limiting[i][x]))
            conc.append([vals[i], vals[i], conc0[0] - vals[i], conc0[0] - vals[i]])
        temp = 25         #Specific to this reaction, reaction occurs at room temperature
        pressure =  1     #Specific to this reaction

        for i in range(len(simulation_state.kinetics)): #conc is the list of lists of concentrations of reactor species. its length is the number of reactors.
            simulation_state.kinetics[i].temperature = temp
            simulation_state.kinetics[i].pressure = pressure
            while(len(simulation_state.kinetics[i].mole_fraction) < len(conc[i])):
                simulation_state.kinetics[i].mole_fraction.append(float(0))
            for j in range(len(conc[i])):
                simulation_state.kinetics[i].mole_fraction[j] = (conc[i][j])
        return simulation_state



