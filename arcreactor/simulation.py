import numpy as np
import datetime as dt
import time
#import scipy.integrate as si
from protobufs.kinetics_pb2 import *
import math

'''
We consider "Hydrolysis of ester" - A psuedo first order chemical reaction.
All reactors are equally sized and operate at the same temperature.
C2H5COOCH3 + H2O --> CH3COOH + C2H5OH
'''

class Simulation:
    '''Controls simulation of objects'''
    def __init__(self, start_time):
        #For demo purposes, the values are fixed
        self.reactor_number = 0
        self.reactor_volumes = np.asarray([2, 2, 2])            # dm3
        self.volumetric_feed_rates = np.asarray([0.5, 0.5])     # dm3/s
        self.molar_feed_rate = np.asarray([1, 1])               # mol/s
        self.start_time = start_time
        self.time = 0

    def edge_list(self, vision_state):
        ''' Reads labels from the vision protobuf and makes a dictionary edge directions showing reactor connections'''
        edge_list = {}

        for j in range(len(vision_state.edges)):
            if(edge_list[vision_state.edges[j]]):
                edge_list[vision_state.edges[j].idA].append(vision_state.edges[j].idB)
            else:
                edge_list[vision_state.edges[j].idA] = [vision_state.edges[j].idB]
        
        return edge_list


    def calculate(self, simulation_state, graph):
        '''Does the actual simulation for number of objects specified by the protobuf '''

        edge_list = self.edge_list(graph)

        R = 8.314      # Universal gas constant (kJ/kmol K)
        T = 298        # Temperature (K), Specific to this reaction, reaction occurs at room temperature
        P =  1         # Pressure (atm), Specific to this reaction
        
        if(len(simulation_state.kinetics) == 0):
            return simulation_state
        if(self.reactor_number != len(simulation_state.kinetics)):  #reset simulation when no of reactors change
            self.reactor_number = len(simulation_state.kinetics)
            self.start_time = simulation_state.time
        conc = []
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        tau = self.reactor_volumes[0] / self.volumetric_feed_rates[0]
        
        k = math.exp(-E_act / (R * T))      # Rate constant
        
        self.time = simulation_state.time - self.start_time
        conc_limiting, vals = [conc0[0]], []
        
        for i in range(len(simulation_state.kinetics)):
            if(simulation_state.kinetics[i].label == 'CSTR'):
                conc_limiting.append(self.cstr(conc_limiting[i-1], tau = tau, k = k))
            elif (simulation_state.kinetics[i].label == 'PFR'):
                conc_limiting.append(self.pfr(conc_limiting[i-1], tau = tau, k = k))

            vals.append(conc_limiting[i])
            conc.append([vals[i], vals[i], conc0[0] - vals[i], conc0[0] - vals[i]])
        
        #conc is the list of lists of concentrations of reactor species. its length is the number of reactors.
            simulation_state.kinetics[i].temperature = T
            simulation_state.kinetics[i].pressure = P
            while(len(simulation_state.kinetics[i].mole_fraction) < len(conc[i])):
                simulation_state.kinetics[i].mole_fraction.append(float(0))
            for j in range(len(conc[i])):
                simulation_state.kinetics[i].mole_fraction[j] = (conc[i][j])
        return simulation_state

    def cstr(self, initial_conc, tau = 1, k = 0.1):
        '''Calculates concentrations for a first order reaction in a CSTR.
        
        Parameters
        ----------
        initial_conc : float
                    Concentration of limiting reactant entering the reactor
        tau : float
            Denotes residence time of the reactor
        k : float
            Denotes the reaction constant for the reaction
        
        Returns
        -------
        float
            Final concentration of the limiting reactor when it leaves the reactor
        '''
        out_conc_lr = initial_conc / ((1 + tau * k))
        return out_conc_lr

    def pfr(self, initial_conc, tau = 1, k = 0.1):
        '''Calculates concentrations for a first order reaction in a PFR.

        Parameters
        ----------
        initial_conc : float
                    Concentration of limiting reactant entering the reactor
        tau : float
            Denotes residence time of the reactor
        k : float
            Denotes the reaction constant for the reaction

        Returns
        -------
        float
                Final concentration of the limiting reactor when it leaves the reactor
        '''
        out_conc_lr = initial_conc * (1 - math.exp**(-k * tau))     
        return out_conc_lr
