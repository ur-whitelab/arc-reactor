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
        ''' Reads labels from the vision protobuf and makes two dictionaries which record inward and outward connections respectively of reactors'''
        edge_list_in = {}
        edge_list_out = {}
        
        for j in range(len(vision_state.edges)):
            if(vision_state.edges[j].idA not in edge_list_out.keys()):
                edge_list_out[vision_state.edges[j].idA] = []
            if(vision_state.edges[j].idB not in edge_list_out.keys()):
                edge_list_out[vision_state.edges[j].idB] = []
            if(vision_state.edges[j].idA not in edge_list_in.keys()):
                edge_list_in[vision_state.edges[j].idA] = []
            if(vision_state.edges[j].idB not in edge_list_in.keys()):
                edge_list_in[vision_state.edges[j].idB] = []
            edge_list_out[vision_state.edges[j].idA].append(vision_state.edges[j].idB)
            edge_list_in[vision_state.edges[j].idB].append(vision_state.edges[j].idA)
        
        return edge_list_in, edge_list_out


    def calculate(self, simulation_state, graph):
        '''Does the actual simulation for number of objects specified by the protobuf '''

        edge_list_in, edge_list_out = self.edge_list(graph)

        if(len(simulation_state.kinetics) == 0):
            return simulation_state
        if(self.reactor_number != len(simulation_state.kinetics)):  #reset simulation when no of reactors change
            self.reactor_number = len(simulation_state.kinetics)
            self.start_time = simulation_state.time
        
        R = 8.314      # Universal gas constant (kJ/kmol K)
        T = 298        # Temperature (K), Specific to this reaction, reaction occurs at room temperature
        P =  1         # Pressure (atm), Specific to this reaction
        
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        tau = self.reactor_volumes[0] / self.volumetric_feed_rates[0]      # seconds
        e_act = 100     #kJ/kmol, yet to use a model to calc e_act as per reaction
        k = math.exp(-e_act / (R * T))      # Rate constant, time dependence needs to be added
        
        self.time = simulation_state.time - self.start_time
        conc_in, reactor_type, conc_limiting = {}, {}, {}

        for i in range(len(simulation_state.kinetics)):
            conc_in[simulation_state.kinetics[i].id] = [conc0[0]]
            reactor_type[simulation_state.kinetics[i].id] = simulation_state.kinetics[i].label
            conc_limiting[simulation_state.kinetics[i].id] = []

        for i in edge_list_out:
            if(edge_list_in[i] == []):
                conc_limiting[i] = calc_conc(conc0, reactor_type[i])
            elif(len(edge_list_in[i]) == 1):
                conc_limiting[i] = calc_conc(conc_in[edge_list_in[i]], reactor_type[i])
            else:
                sum_conc, vol = 0, 0
                for j in edge_list_in[i]:
                    sum_conc += conc_in[j] * self.volumetric_feed_rates[i]
                    vol += self.volumetric_feed_rates[i]
                conc_limiting[i] = calc_conc(sum_conc/vol, reactor_type[i])

            conc = [conc_limiting[i], conc_limiting[i], conc0[0] - conc_limiting[i], conc0[0] - conc_limiting[i]]
            conc_in[i] = conc_limiting[i]
            #conc is the list of lists of concentrations of chemical species. It's length is the number of reactors.
            simulation_state.kinetics[i].temperature = T
            simulation_state.kinetics[i].pressure = P
            while(len(simulation_state.kinetics[i].mole_fraction) < len(conc[i])):
                simulation_state.kinetics[i].mole_fraction.append(float(0))
            for j in range(len(conc[i])):
                simulation_state.kinetics[i].mole_fraction[j] = conc
        return simulation_state

                
        def calc_conc(initial_conc, reactor_type):
            if(reactor_type == 'CSTR'):
                conc_limiting = self.cstr(initial_conc, tau = tau, k = k)
            elif(reactor_type == 'PFR'):
                conc_limiting = self.pfr(initial_conc, tau = tau, k = k)
            return conc_limiting


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
                Final concentration of the limiting reactant when it leaves the reactor
        '''
        out_conc_lr = initial_conc * (math.exp**(-k * tau))
        return out_conc_lr
