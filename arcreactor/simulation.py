import numpy as np
import datetime as dt
import time
import scipy.integrate as si
from .protobufs.kinetics_pb2 import *
import math
import sys

'''
We consider a pseudo first order reversible chemical reaction which is equilibrium limited.
All reactors are equally sized and participants decide the temperature at which they are operated.
A + B <--> C + D

Update: Now using Sabatier equation.
CO2 + 4H2 -> CH4 + 2H2O
'''

class Simulation:
    '''Controls simulation of objects'''
    def __init__(self, start_time):
        #For demo purposes, the values are fixed
        self.reactor_number = 0
        self.reactor_volume = 200            # m3
        self.volumetric_feed_rates = np.array([10, 10])     # m3/s
        self.molar_feed_rate = np.array([1, 1])           # mol/s
        self.start_time = start_time
        self.graph_time = 0
        self.time = 0
        self.edge_list_in = {}
        self.edge_list_out = {}
        self.vol_out_rates = {}
        self.vol_in_rates = {}
        self.connected_to_source = False
        self.edge_list_changed = False
        self.conc0 = self.molar_feed_rate[0]/self.volumetric_feed_rates[0]  # mol/dm3
        self.start_plotting = False #flag to start the plots
        self.restart_plots = False
        #stoichiometry
        self.a = 1
        self.b = 4
        self.c = 1
        self.d = 2
        self.ready_flags = {}#these are for tracking when PFRs are finished reacting
        self.ready_flags[0] = True #Source is always ready!
        self.done_times = {}#A reactor only starts outputting if all its incoming edges are done
        self.done_times[0] = 0.0 #these are times so they should be floats


    def update_edge_list(self, graph):
        ''' Reads labels from the vision protobuf and makes a dictionary which records inward connections of each reactor'''
        self.edge_list_out[0] = []
        for key in graph.nodes:
            node = graph.nodes[key]
            if((node.id not in self.edge_list_in and not node.delete) and (node.id != 999 and node.id != 0)):#don't add for the conditions or source nodes; they never take in
                self.edge_list_in[node.id] = []#new ID, make new lists for it
                self.edge_list_changed = True
            if((node.id not in self.edge_list_out and not node.delete) and (node.id != 999 and node.id != 0)):#don't add for the conditions or source nodes; they never take in
                self.edge_list_out[node.id] = []
                self.edge_list_changed = True
            elif(node.delete):
                self.edge_list_changed = True
                for edgekey in self.edge_list_in:
                    if(node.id in self.edge_list_in[edgekey]):
                        self.edge_list_in[edgekey].remove(node.id)
                    if(edgekey == node.id):
                        if(0 in self.edge_list_in[node.id]):
                            self.connected_to_source = False
                        self.edge_list_in[edgekey] = []#empty it
                for edgekey in self.edge_list_out:
                    if(node.id in self.edge_list_out[edgekey]):
                        self.edge_list_out[edgekey].remove(node.id)
                    if(edgekey == node.id):
                        if(0 in self.edge_list_out[node.id]):
                            self.connected_to_source = False
                        self.edge_list_out[edgekey] = []#empty it

        for key in graph.edges:
            edge = graph.edges[key]
            if ((edge.idB in self.edge_list_in) and (edge.idA not in self.edge_list_in[edge.idB]) or len(self.edge_list_in[edge.idB]) == 0):#append if it's a new node to this one
                self.edge_list_in[edge.idB].append(edge.idA)
                self.edge_list_changed = True
            if ((edge.idA in self.edge_list_out) and (edge.idB not in self.edge_list_out[edge.idA]) or len(self.edge_list_out[edge.idA]) == 0):#append if it's a new node from this one
                self.edge_list_out[edge.idA].append(edge.idB)
                self.edge_list_changed = True
            if(edge.idA == 0):#source
                self.connected_to_source = True
                self.edge_list_out[0].append(edge.idB)

    def update_out_rates(self, id):
        '''Called recursively to calculate volumetric out rates. Not displayed.'''
        if(id == 0):
            self.vol_out_rates[id] = self.volumetric_feed_rates[0]
        else:
            vol_in_sum = 0.0
            for node in self.edge_list_in[id]:
                if(node in self.vol_out_rates):
                    val = self.vol_out_rates[node]
                else:
                    val=0.0
                vol_in_sum += val
            self.vol_in_rates[id] = vol_in_sum
            self.vol_out_rates[id] = vol_in_sum / max(len(self.edge_list_out[id]), 1.0)
        if(len(self.edge_list_out[id]) == 0):
            return
        for key in self.edge_list_out[id]:
            self.update_out_rates(key)


        #print('self.edge_list_in is {}\n'.format(self.edge_list_in))



    def add_delete_protobuf_objects(self, simulation_state, graph):
        '''Add and delete kinetics objects from the kinetics protobuf, and assign id, label, and temperature to new objects'''
        #delete the whole list each time and create a new one
        #TODO: need to fix this so we can re-add reactors and connect them again. #UPDATE: Reactors can be re-added but are not reset unless source connection is removed.
        length = len(simulation_state.kinetics)
        for i in range(length):
            del simulation_state.kinetics[-1]
        for key in graph.nodes:
            node = graph.nodes[key]
            if (node.id != 999 and node.id != 0):
                if(node.delete is not True):
                    simulation_state.kinetics.add() #always just append
                    simulation_state.kinetics[-1].label = node.label #always get the current last kinetics object
                    simulation_state.kinetics[-1].id = node.id
                    self.ready_flags[node.id] = False #start not ready
                    if (len(node.weight) > 0):
                        simulation_state.kinetics[-1].temperature = node.weight[0]#T is always the first in this repeat field
                    else:
                        simulation_state.kinetics[-1].temperature = 400  #default
        return simulation_state


    def calc_conc(self, initial_conc, reactor_type, k_eq, k, id):
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates #molar_feed_rate is a list
        if(reactor_type == 'cstr'):
            conc_limiting, ready = self.cstr(initial_conc, t=self.time, k_eq = k_eq, k = k)
        elif(reactor_type == 'pfr'):
            conc_limiting, ready = self.pfr(initial_conc, t=self.time, k_eq = k_eq, k = k, done_time = self.done_times[id])
        else:
            conc_limiting = conc0[0]
        #print(conc_limiting)
        return (conc_limiting, ready)

    def calc_outputs(self, id, simulation_state, R, P):
        '''RECURSIVELY calculate output concentrations for each node in the graph. ALWAYS call with id 0 first!'''
        found = False #assume it's not there
        for kinetics in simulation_state.kinetics:#look for node with incoming connection from this ID
            if(id in self.edge_list_in[kinetics.id]):
                found = True
                #found it! set output concentrations and recurse
                if(kinetics.temperature != 0):
                    T = kinetics.temperature
                    e_act = 47000  #kJ/kmol
                    k_eq = 100000 * math.e ** (-33.78*(T-298)/T)    #equilibrium constant
                    k = 5*10**6 * math.exp(-e_act / (R * T))      # Rate constant, time dependence needs to be added
                    #find the limiting concentration for the ith reactor
                    #conc_limiting = self.calc_conc(sum([conc_out[idx] for idx in self.edge_list_in[i]]), kinetics.label, kinetics.id, k_eq, k)
                    conc_out_sum = 0.0
                    conc_product = 0.0
                    vol_out_sum = 0.0
                    all_incoming_ready = True
                    for idx in self.edge_list_in[kinetics.id]:
                        if(self.ready_flags[idx] == False):
                            all_incoming_ready = False
                    if(all_incoming_ready):
                        max_done_time_in = 0.0
                        for idx in self.edge_list_in[kinetics.id]:
                            conc_out_sum += self.conc_out_reactant[idx] * self.vol_out_rates[kinetics.id] #len(self.edge_list_in[kinetics.id])
                            #keeping track of product coming out of previous reactors
                            conc_product += self.conc_out_product[idx] * self.vol_out_rates[kinetics.id]
                            vol_out_sum += self.vol_out_rates[kinetics.id]
                            max_done_time_in = max(max_done_time_in, self.done_times[idx])
                        conc_out_sum /= vol_out_sum# (C1V1 + C2V2)/(V1+V2) = C_final
                        conc_product /= vol_out_sum
                        if(kinetics.label == 'cstr'):
                            self.done_times[kinetics.id] = 0.0
                        elif(kinetics.label == 'pfr'):
                            self.done_times[kinetics.id] = max_done_time_in + self.reactor_volume/self.vol_in_rates[kinetics.id]
                        #print('conc_out_sum is {} and conc_product is {} for id {}'.format(conc_out_sum, conc_product, kinetics.id))
                        conc_limiting, self.ready_flags[kinetics.id] = self.calc_conc(initial_conc = conc_out_sum, reactor_type=kinetics.label,  k_eq=k_eq, k=k, id=kinetics.id)
                        #print('Conc limiting for reactor id {} is {}'.format(kinetics.id, conc_limiting))

                        #conc = [conc_limiting, conc_limiting, conc_product + (conc_out_sum - conc_limiting), conc_product + (conc_out_sum - conc_limiting)]
                        #print('conc is {}'.format(conc))
                        self.conc_out_reactant[kinetics.id] = conc_limiting
                        self.conc_out_product[kinetics.id] = conc_product + conc_out_sum - conc_limiting   #taking into account the existing conc of products
                    else:#Do NOT output until ready
                        self.conc_out_reactant[kinetics.id] = 0.0
                        self.conc_out_product[kinetics.id] = 0.0
                        self.ready_flags[kinetics.id] = False
                    self.calc_outputs(kinetics.id, simulation_state, R, P)#now that this node has its outputs set, go on to the ones it outputs to

        if(not found):
            return


    async def calculate(self, simulation_state, graph):
        '''The actual simulation for number of objects specified by the protobuf '''
        #graph = graph # update the graph object when we get it (see controller.py)
        self.update_edge_list(graph)
        self.update_out_rates(0)#ONLY call this after update_edge_list() is done, and ONLY with id == 0
        if(self.graph_time % 100 ==0):
            print('update_out_rates was called and now out rates are: {}'.format(self.vol_out_rates))
        simulation_state = self.add_delete_protobuf_objects(simulation_state, graph)
        if(not self.connected_to_source ):
            self.start_plotting = False
            if( self.graph_time % 20 == 0):
                print('NOT CONNECTED TO SOURCE')
            return simulation_state


        if(self.reactor_number != len(simulation_state.kinetics)):
            self.edge_list_changed = True

        if(self.edge_list_changed):  #reset simulation when edges change
            self.reactor_number = len(simulation_state.kinetics)
            self.start_time = self.graph_time
            self.edge_list_changed = False
            self.restart_plots = True
            self.start_plotting = True
            #print('reactors = {}'.format(self.reactor_number))

        R = 8.314      # Universal gas constant (kJ/kmol K)
        P = 1          # Pressure (atm)


        self.time = self.graph_time - self.start_time
        #print('Time is {}'.format(self.time))
        self.conc_out_reactant, self.conc_out_product = {0:self.conc0}, {0:0}

        for kinetics in simulation_state.kinetics:
            self.conc_out_reactant[kinetics.id] = 0
            self.conc_out_product[kinetics.id] = 0
            kinetics.pressure = P
        self.calc_outputs(id = 0, simulation_state = simulation_state, R = R, P = P)
        #start with ID zero to look for nodes connected to source

        #print('simulation_state.kinetics is {}'.format(simulation_state.kinetics))
        #print('the keys for conc_out are {}, and the keys for self.edge_list_in are {}'.format(conc_out.keys(), self.edge_list_in.keys()))
        #sys.stdout.flush()
        for kinetics in simulation_state.kinetics:#TODO: change this loop to work "smartly" from source to ends... Recursion?
            i = kinetics.id
            conc_limiting = self.conc_out_reactant[kinetics.id]
            conc_out_product = self.conc_out_product[kinetics.id]  #taking into account the existing conc of products
            conc = [conc_limiting, self.b / self.a * conc_limiting, conc_out_product, self.d / self.c * conc_out_product]
            for j in range(len(conc)):
                kinetics.mole_fraction.append(float(conc[j]))
                #if(simulation_state.time %5 == 0):
                    #print('The {}th mole fractions are {}'.format(i, kinetics.mole_fraction))
        simulation_state.time = self.time
        return simulation_state


    def cstr(self, initial_conc, t, k_eq = 5, k = 0.1):
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
           return -k* initial_conc * (1 - (1 + self.c / self.a /k_eq)*conv)
        rv = self.reactor_volume #m3
        fa0 = self.molar_feed_rate[0] #mol/s
        conversion = min(rv * k * initial_conc/(fa0 + k* rv * initial_conc + (self.c / self.a * k * rv * initial_conc)/k_eq), 1.0)
        #print('Conversion from cstr is {} at {}'.format(conversion, self.time))
        out_conc_lr = initial_conc*(1.0 - conversion)
        ready = True #CSTR is instantaneous
        #print('A left in cstr is {}'.format(out_conc_lr))
        return (out_conc_lr, ready)

    def pfr(self, initial_conc, t, k_eq = 5, k = 0.1, done_time = None):
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
        if(done_time is None):
            done_time = self.reactor_volume/self.volumetric_feed_rates[1]
        def rate(conv, t):#must pass in start time.
            return k * initial_conc*(1 - (1 + self.c / self.a /k_eq) * conv)
        #conversion = min(si.odeint(rate, 0.0001, np.arange(0, 3600, 3600*25)), 1.0)  # ~25fps
        time = min(done_time, self.time/3.25)
        ready = False
        if(self.time/3.25 >= done_time):
            ready = True
        conversion = min((k_eq - k_eq * math.exp( -time * k * (self.c/self.a + k_eq) / k_eq))/(k_eq + self.c/self.a), 1) #using derived formula; divide by 3.25 to go from FPS to accelerated seconds
        #print('Conversion from pfr is {} at {}'.format(conversion, t))
        out_conc_lr = initial_conc*(1.0 - conversion)
        #print('A left in pfr is {}'.format(out_conc_lr))
        return (out_conc_lr, ready)
