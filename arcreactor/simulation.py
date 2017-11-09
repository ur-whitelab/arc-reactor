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
        self.volumetric_feed_rates = np.array([10, 10])     # m3/s
        self.molar_feed_rate = np.array([1, 1])           # mol/s
        self.start_time = start_time
        self.time = 0

    def edge_list(self, graph):
        ''' Reads labels from the vision protobuf and makes two dictionaries which record inward and outward connections respectively of reactors'''
        edge_list_in = {}
        edge_list_out = {}

        for i in range(len(graph.nodes)):
            edge_list_in[graph.nodes[i].id] = []
            edge_list_out[graph.nodes[i].id] = []

        if(len(graph.edges) > 0):
            for i in range(len(graph.edges)):
                #should already be allocated?
                if graph.edges[i].idA not in edge_list_out:
                    edge_list_out[graph.edges[i].idA] = []
                #if graph.edges[i].idA not in edge_list_in:
                #    edge_list_in[graph.edges[i].idA] = []
                if graph.edges[i].idB not in edge_list_in:
                    edge_list_in[graph.edges[i].idB] = []
                #if graph.edges[i].idB not in edge_list_out:
                #    edge_list_out[graph.edges[i].idB] = []

                edge_list_out[graph.edges[i].idA].append(graph.edges[i].idB)
                edge_list_in[graph.edges[i].idB].append(graph.edges[i].idA)
        #print(edge_list_in)
        return edge_list_in, edge_list_out


    def add_delete_protobuf_objects(self, simulation_state, graph):
        '''Add and delete kinetics objects from the kinetics protobuf, and assign id, label, and temperature to new objects'''
        T = 423
        #delete the whole list each time and create a new one
        for i in range(len(simulation_state.kinetics)):
            del simulation_state.kinetics[-1]
        for i in range(len(graph.nodes)):
            if(graph.nodes[i].delete):
                simulation_state.kinetics.add() #empty placeholders
            elif(graph.nodes[i].id != 999): #non-source non-empty non-conditions nodes
                simulation_state.kinetics.add()
                simulation_state.kinetics[i].label = graph.nodes[i].label
                print('the label for this node is {} and its id is {}'.format(graph.nodes[i].label, graph.nodes[i].id))#this ID is coming out as 0. Why?
                simulation_state.kinetics[i].id = graph.nodes[i].id
                simulation_state.kinetics[i].temperature = T  #default
        simulation_state.kinetics.add()
        simulation_state.kinetics[0].id = 0
        simulation_state.kinetics[0].label = 'source' #Unity is seeing a kinetics object with ID 0 for some reason...
        return simulation_state


    def calc_conc(self, initial_conc, reactor_type, i, k_eq, k):
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates
        if(reactor_type == 'cstr'):
            conc_limiting = self.cstr(initial_conc, self.time, i, k_eq = k_eq, k = k)
        elif(reactor_type == 'pfr'):
            conc_limiting = self.pfr(initial_conc, self.time, k_eq = k_eq, k = k)
        else:
            conc_limiting = conc0[0]
        #print(conc_limiting)
        return conc_limiting

    def calculate(self, simulation_state, graph):
        '''The actual simulation for number of objects specified by the protobuf '''
        #graph = graph # update the graph object when we get it (see controller.py)
        if(len(graph.edges) == 0 or len(graph.nodes) == 0): #check if there are any nodes and edges
            return simulation_state
        connected_to_source = False
        for i in range(len(graph.edges)):   #start simulation only if one of the reactors is connected to the source
            if(graph.edges[i].labelA == 'source'):
                print('CONNECTED TO SOURCE')
                connected_to_source = True
        if(not connected_to_source):
            print('NOT CONNECTED TO SOURCE')
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
        conc_out, reactor_type, conc_limiting = {0:conc0[0]}, {}, {}

        for kinetics in simulation_state.kinetics:
            conc_out[kinetics.id] = conc0[0]
            #record concentration coming out of the reactors
            reactor_type[kinetics.id] = kinetics.label

        print('simulation_state.kinetics is {}'.format(simulation_state.kinetics))
        print('the keys for conc_out are {}, and the keys for edge_list_in are {}, and for edge_list_out they are {}'.format(conc_out.keys(), edge_list_in.keys(), edge_list_out.keys()))
        for kinetics in simulation_state.kinetics:
            i = kinetics.id
            if(kinetics.temperature != 0):
                T = kinetics.temperature
                e_act = 25  #kJ/kmol
                k_eq = 100000 * math.e ** (-33.78*(T-298)/T)    #equilibrium constant
                k = math.exp(-e_act / (R * T))      # Rate constant, time dependence needs to be added
                #find the limiting concentration for the ith reactor
                conc_limiting = self.calc_conc(sum([conc_out[idx] for idx in edge_list_in[i]]), kinetics.label, kinetics.id, k_eq, k)
                #print('Conc limiting is {}'.format(conc_limiting))

                conc = [conc_limiting, conc_limiting, (conc0[0] - conc_limiting), (conc0[0] - conc_limiting)]
                #print('conc is {}'.format(conc))
                conc_out[kinetics.id] = conc_limiting
                #conc is the list of lists of concentrations of chemical species. It's length is the number of reactors.
                kinetics.temperature = T
                kinetics.pressure = P
                j=0
                while(len(kinetics.mole_fraction) < len(conc)):
                    kinetics.mole_fraction.append(float(0))
                for j in range(len(conc)):
                    kinetics.mole_fraction[j] = float(conc[j])
                if(simulation_state.time %5 == 0):
                    print('The {}th kinetics is {}'.format(i, kinetics))

        return simulation_state


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
        #def rate(conv, t):
        #   return -k* initial_conc * (1 - (1 + 1/k_eq)*conv)

        #if(i==1):       #CSTR is a steady state reactor, but over a long time we'd get a constant value
        #    conv = si.odeint(rate, 0.0001, np.arange(0, 3600, 3600*25))
        #    out_conc_lr = initial_conc * (1 - conv[int(t)])
        #else:
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
