import numpy as np
import datetime as dt
import time
import scipy.integrate as si
from protobufs.kinetics_pb2 import *
import math

'''
We consider a psuedo first order reversible chemical reaction which is equilibrium limited.
All reactors are equally sized and participants decide the temperature at which they are operated.
A + B <--> C + D
'''

class Simulation:
    '''Controls simulation of objects'''
    def __init__(self, start_time):
        #For demo purposes, the values are fixed
        self.reactor_number = 0
        self.reactor_volume = 20            # m3
        self.volumetric_feed_rates = np.asarray([10, 10])     # m3/s
        self.molar_feed_rate = np.asarray([1, 1])           # mol/s
        self.start_time = start_time
        self.time = 0

    def edge_list(self, vision_state):
        ''' Reads labels from the vision protobuf and makes two dictionaries which record inward and outward connections respectively of reactors'''
        edge_list_in = {}
        edge_list_out = {}
        
        for i in range(len(vision_state.nodes)):
            edge_list_in[vision_state.nodes[i].id] = []
            edge_list_out[vision_state.nodes[i].id] = []

        if(len(vision_state.edges) > 0):
            for j in range(len(vision_state.edges)):
                edge_list_out[vision_state.edges[j].idA].append(vision_state.edges[j].idB)
                edge_list_in[vision_state.edges[j].idB].append(vision_state.edges[j].idA)
        
        return edge_list_in, edge_list_out


    def add_delete_protobuf_objects(self, simulation_state, graph):
        '''Add and delete kinetics objects from the kinetics protobuf, and assign id, label, and temperature to new objects'''
        j = 0       #for keeping a count of simulation_state.kinetics objects
        T = 423
        for i in range(len(graph.nodes)):
            if(graph.nodes[i].label == 'source'):
                continue
            if(graph.nodes[i].delete is True):
                continue
            else:
                if(simulation_state.kinetics[j] is True):
                    simulation_state.kinetics[j].label = graph.nodes[i].label
                    simulation_state.kinetics[j].id = graph.nodes[i].id
                    if(graph.nodes[i].weight == 0):
                        simulation_state.kinetics[j].temperature = T  #default
                    else:
                        simulation_state.kinetics[j].temperature = graph.nodes[i].weight
                    j += 1
                else:
                    simulation_state.kinetics.add()
                    simulation_state.kinetics[j].label = graph.nodes[i].label
                    simulation_state.kinetics[j].id = graph.nodes[i].id
                    if(graph.nodes[i].weight == 0):
                        simulation_state.kinetics[j].temperature = T  #default
                    else:
                        simulation_state.kinetics[j].temperature = graph.nodes[i].weight
                    j += 1
        if(len(simulation_state.kinetics) > j): #delete extra kinetics objects
            for i in range(len(simulation_state.kinetics) - j):
                del simulation_state.kinetics[-1]
        return simulation_state


    def calculate(self, simulation_state, graph):
        '''The actual simulation for number of objects specified by the protobuf '''

        if(len(graph.edges) == 0 or len(graph.nodes) == 0): #check if there are any nodes and edges
            return simulation_state
        for i in range(len(graph.edges)):   #start simulation only if one of the reactors is connected to the source
            if(graph.edges[i].labelA == 'source'):
                break
            else:
                return simulation_state

        simulation_state = self.add_delete_protobuf_objects(simulation_state, graph)
        edge_list_in, edge_list_out = self.edge_list(graph) 

        if(self.reactor_number != len(simulation_state.kinetics)):  #reset simulation when no of reactors change
            self.reactor_number = len(simulation_state.kinetics)
            self.start_time = simulation_state.time
        
        R = 8.314      # Universal gas constant (kJ/kmol K)
        P = 1          # Pressure (atm)
        
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates  # mol/dm3
        #tau = self.reactor_volume / self.volumetric_feed_rates[0]      # seconds

        def find_activation_energy(T):
            '''Interpolate through known values of activation energy and return activation energy corresponding to temperature T'''
            temp = [298, 348, 398, 443, 448, 498, 548, 598, 648, 698, 748, 798, 848, 898, 948, 998]
            e_act = [20.25, 21.87, 24.01, 25, 26.69, 29.96, 33.85, 38.41, 43.66, 49.65, 56.41, 63.98, 72.4, 81.71, 91.94, 103.13]
            #kJ/kmol, yet to use a model to calc e_act as per reaction 
            e_act1 = np.interp(T, temp, e_act)
            return e_act1

        self.time = simulation_state.time - self.start_time
        conc_out, reactor_type, conc_limiting = {}, {}, {}

        for i in range(len(simulation_state.kinetics)):
            if(simulation_state.kinetics[i].label == 'source'):
                conc_out[simulation_state.kinetics[i].id] = [conc0[0]]
            else:    
                conc_out[simulation_state.kinetics[i].id] = []      #record concentration coming out of the reactors
            reactor_type[simulation_state.kinetics[i].id] = simulation_state.kinetics[i].label
            conc_limiting[simulation_state.kinetics[i].id] = []     

        for i in range(len(simulation_state.kinetics)):
            T = simulation_state.kinetics[i].temperature
            e_act = find_activation_energy(T)
            k_eq = 100000 * math.e ** (-33.78*(T-298)/T)    #equilibrium constant
            k = math.exp(-e_act / (R * T))      # Rate constant, time dependence needs to be added
            if(len(edge_list_in[i]) == 1):
                conc_limiting[i] = calc_conc(conc_out[edge_list_in[i]], reactor_type[i], simulation_state.kinetics[i].id, k_eq, k)
            else:
                sum_conc, vol = 0, 0
                for j in edge_list_in[i]:
                    sum_conc += conc_out[j] * self.volumetric_feed_rates[i]
                    vol += self.volumetric_feed_rates[i]
                conc_limiting[i] = calc_conc(sum_conc/vol, reactor_type[i], simulation_state.kinetics[i].id, k_eq, k)

            conc = [conc_limiting[i], conc_limiting[i], conc0[0] - conc_limiting[i], conc0[0] - conc_limiting[i]]
            conc_out[i] = conc_limiting[i]
            #conc is the list of lists of concentrations of chemical species. It's length is the number of reactors.
            simulation_state.kinetics[i].temperature = T
            simulation_state.kinetics[i].pressure = P
            while(len(simulation_state.kinetics[i].mole_fraction) < len(conc[i])):
                simulation_state.kinetics[i].mole_fraction.append(float(0))
            for j in range(len(conc)):
                simulation_state.kinetics[i].mole_fraction[j] = conc[j]
        return simulation_state

        def calc_conc(initial_conc, reactor_type, i, k_eq, k):
            if(reactor_type == 'CSTR'):
                conc_limiting = self.cstr(initial_conc, self.time, i, k_eq = k_eq, k = k)
            elif(reactor_type == 'PFR'):
                conc_limiting = self.pfr(initial_conc, self.time, k_eq = k_eq, k = k)
            return conc_limiting


    def cstr(self, initial_conc, i, t, k_eq = 5, k = 0.1):
        '''Calculates concentrations for a first order reaction in a CSTR.
        
        Parameters
        ----------
        initial_conc : float
                    Concentration of limiting reactant entering the reactor
        i : int
            Reactor id
        k_eq : float
            Denotes equilibrium concentration of the reaction
        k : float
            Denotes the reaction constant for the reaction
        
        Returns
        -------
        float
            Final concentration of the limiting reactor when it leaves the reactor
        '''
        def rate(conv, t):
            return -k* initial_conc * (1 - (1 + 1/k_eq)*conv)

        if(i==1):       #CSTR is a steady state reactor, but over a long time we'd get a constant value
            conv = si.odeint(rate, 0.0001, np.arange(0, 3600, 3600*25))
            out_conc_lr = initial_conc * (1 - conv[int(t)])
        else:
            rv = 20 #m3
            fa0 = 1 #mol/s
            v0 = 10 #m3/s
            conversion = rv*k*initial_conc/(fa0+k*rv*initial_conc+(k*rv*initial_conc)/k_eq)
            out_conc_lr = initial_conc*(1 - conversion)
        return out_conc_lr

    def pfr(self, initial_conc, t, k_eq = 5, k = 0.1):
        '''Calculates concentrations for a first order reaction in a PFR.

        Parameters
        ----------
        initial_conc : float
                    Concentration of limiting reactant entering the reactor
        t: int
            Time in simulation at which concentration needs to be calculated
        k_eq : float
            Denotes equilibrium concentration of the reaction
        k : float
            Denotes the reaction constant for the reaction

        Returns
        -------
        float
                Final concentration of the limiting reactant when it leaves the reactor
        '''
        def rate(conv, t):
            return k * initial_conc*(1 - (1+ 1/k_eq) * conv)

        rv = 20 #m3
        fa0 = 1 #mol/s
        v0 = 10 #m3/s
        conversion = si.odeint(rate, 0.0001, np.arange(0, 3600, 3600*25))  # ~25fps
        out_conc_lr = initial_conc*(1 - conversion[int(t)])
        return out_conc_lr
