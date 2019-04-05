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
Sabatier is a gas phase reaction and hence, cannot be conducted in a CSTR
So now using: 2 HCl(aq) + Mg(OH)2(aq) <--> 2 H2O(ℓ) + MgCl2(aq)
This reaction is aqueous and not exciting. So, switching now to friedel-craft reaction.
Benzene + 3C2H5Br <--> 1,3,5-triethylbenzene + 3 HBr - First order in C2H5Br
'''

class Simulation:
    '''Controls simulation of objects'''
    def __init__(self, start_time):
        #For demo purposes, the values are fixed
        #self.chemical_species = ['Benzene', 'EtBr', 'TEB', 'HBr']
        self.chemical_species = ['A', 'B']#, '', '']
        self.reactor_number = 0
        self.volumetric_feed_rates = np.array([1., 1.])   # dm3/s, i.e. L/s
        self.molar_feed_rate = np.array([1., 1.])           # mol/s
        self.start_time = start_time
        self.graph_time = 0
        self.time = 0
        self.edge_list_in = {}
        self.edge_list_out = {0:[]}
        self.vol_out_rates = {0:self.volumetric_feed_rates[0]}
        self.vol_in_rates = {}
        self.connected_to_source = False
        self.edge_list_changed = False
        self.conc0 = self.molar_feed_rate[0]/self.volumetric_feed_rates[0]  # mol/dm3 i.e. mol/L
        self.start_plotting = False #flag to start the plots
        self.restart_plots = False
        #stoichiometry
        self.a = 1 #reactant 1
        self.b = 0 #reactant 2
        self.c = 1 #product 1
        self.d = 0 #product 2
        self.ready_flags = {}#these are for tracking when PFRs are finished reacting
        self.ready_flags[0] = True #Source is always ready!
        self.done_times = {}#A reactor only starts outputting if all its incoming edges are done
        self.done_times[0] = 0.0 #these are times so they should be floats


    def update_edge_list(self, graph):
        ''' Reads labels from the vision protobuf and makes a dictionary which records inward connections of each reactor'''
        for key in graph.nodes:
            node = graph.nodes[key]
            if (node.id not in self.edge_list_in and not node.delete) and (node.id != 999 and node.id != 0):#don't add for the conditions or source nodes; they never take in
                self.edge_list_in[node.id] = []#new ID, make new lists for it
                self.vol_in_rates[node.id] = 0.0
                self.edge_list_changed = True
            if (node.id not in self.edge_list_out and not node.delete) and (node.id != 999):#don't add for the conditions node; it never takes in
                self.edge_list_out[node.id] = []
                self.vol_out_rates[node.id] = 0.0
                self.edge_list_changed = True
            elif node.delete:
                self.edge_list_changed = True
                if node.id in self.vol_out_rates:#if a node is deleted, take it out of the respective dicts
                    self.vol_out_rates.pop(node.id, None)
                if node.id in self.vol_in_rates:
                    self.vol_in_rates.pop(node.id, None)
                for edgekey in self.edge_list_in:
                    if node.id in self.edge_list_in[edgekey]:
                        self.edge_list_in[edgekey].remove(node.id)
                    if edgekey == node.id:
                        if 0 in self.edge_list_in[node.id]:
                            self.connected_to_source = False
                        self.edge_list_in[edgekey] = []#empty it
                for edgekey in self.edge_list_out:
                    if node.id in self.edge_list_out[edgekey]:
                        self.edge_list_out[edgekey].remove(node.id)
                    if edgekey == node.id:
                        self.edge_list_out[edgekey] = []#empty it

        for key in graph.edges:
            edge = graph.edges[key]
            if (edge.idB in self.edge_list_in) and (edge.idA not in self.edge_list_in[edge.idB]) or len(self.edge_list_in[edge.idB]) == 0:#append if it's a new node to this one
                self.edge_list_in[edge.idB].append(edge.idA)
                self.edge_list_changed = True
            if (edge.idA in self.edge_list_out) and (edge.idB not in self.edge_list_out[edge.idA]) or len(self.edge_list_out[edge.idA]) == 0:#append if it's a new node from this one
                self.edge_list_out[edge.idA].append(edge.idB)
                self.edge_list_changed = True
            if edge.idA == 0:#source
                self.connected_to_source = True

    def update_out_rates(self, id):
        '''Called recursively to calculate volumetric out rates. Not displayed.'''
        if(id == 0):
            self.vol_out_rates[id] = self.volumetric_feed_rates[0] / float(len(self.edge_list_out[id]))
        else:
            vol_in_sum = 0.0
            for node in self.edge_list_in[id]:
                if(node in self.vol_out_rates):
                    val = self.vol_out_rates[node]
                else:
                    val=0.0
                vol_in_sum += val
            self.vol_in_rates[id] = vol_in_sum
            self.vol_out_rates[id] = vol_in_sum / max(len(self.edge_list_out[id]), 1)
        if(len(self.edge_list_out[id]) == 0):
            return
        for key in self.edge_list_out[id]:
            self.update_out_rates(key)




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
                        simulation_state.kinetics[-1].pressure = node.weight[1] #actually is volume TODO: change protobuf.
                    else:
                        simulation_state.kinetics[-1].temperature = 393  #default
                        simulation_state.kinetics[-1].pressure = 273 #default; actually is volume. TODO: change protobuf.
        return simulation_state


    def calc_conc(self, initial_conc, initial_conversion, V, reactor_type, k_eq, k, id):
        conc0 = self.molar_feed_rate / self.volumetric_feed_rates #molar_feed_rate is a list
        if(reactor_type == 'cstr'):
            conc_limiting, ready = self.cstr(initial_conc = initial_conc, initial_conversion = initial_conversion,  V = V, k_eq = k_eq, k = k)
        elif(reactor_type == 'pfr'):
            conc_limiting, ready = self.pfr(initial_conc, V = V,  k_eq = k_eq, k = k, done_time = self.done_times[id])
        elif(reactor_type == 'pbr'):
            #conc_limiting, ready = self.pbr(initial_conc, V = V,  k_eq = k_eq, k = k)    for PBR
            conc_limiting, ready = self.pbr(initial_conc, R = 8.314, k_eq = k_eq, k = k)    #for BR
        else:
            conc_limiting = conc0[0]
        return (conc_limiting, ready)

    def calc_outputs(self, id, simulation_state, R):
        '''RECURSIVELY calculate output concentrations for each node in the graph. ALWAYS call with id 0 first!'''
        found = False #assume it's not there
        for kinetics in simulation_state.kinetics:#look for node with incoming connection from this ID
            if(id in self.edge_list_in[kinetics.id]):
                found = True
                #found it! set output concentrations and recurse
                if(kinetics.temperature != 0):
                    T = kinetics.temperature
                    e_act = 47000  #kJ/kmol
                    k_eq = 0.01 * math.exp(15000 /(R * T))#100000 * math.exp(-33.78*(T-298)/T)    #equilibrium constant
                    #k = 5*10**6 * math.exp(-e_act / (R * T))      # Rate constant, time dependence needs to be added
                    k = 100*math.exp(-20000/(R*T))
                    #k_b = #10000*math.exp(-35000/(R*T))
                    V = kinetics.pressure #this is actually volume. TODO: Change protobuf.
                    #find the limiting concentration for the ith reactor
                    #conc_limiting = self.calc_conc(sum([conc_out[idx] for idx in self.edge_list_in[i]]), kinetics.label, kinetics.id, k_eq, k)
                    conc_in_sum = 0.0 # sum of incoming concentrations (conc_out of each incoming rxr)
                    conc_product = 0.0
                    vol_in_sum = 0.0 # sum of incoming concentrations (vol_out_rate of each incoming rxr)
                    all_incoming_ready = True
                    for idx in self.edge_list_in[kinetics.id]:
                        if(self.ready_flags[idx] == False):
                            all_incoming_ready = False
                    if(all_incoming_ready):
                        max_done_time_in = 0.0
                        for idx in self.edge_list_in[kinetics.id]:
                            val = self.vol_out_rates[idx]
                            #concentration of reactants entering the reactor
                            conc_in_sum += self.conc_out_reactant[idx] * val  #len(self.edge_list_in[kinetics.id])
                            #keeping track of product coming out of previous reactors
                            conc_product += self.conc_out_product[idx] * val
                            vol_in_sum += val
                            max_done_time_in = max(max_done_time_in, self.done_times[idx])
                        conc_in_sum /= vol_in_sum# (C1V1 + C2V2)/(V1+V2) = C_final
                        conc_product /= vol_in_sum
                        if(kinetics.label == 'cstr'): #or kinetics.label == 'pbr'):
                            self.done_times[kinetics.id] = max_done_time_in + 0.0
                        elif(kinetics.label == 'pfr'):
                            if self.vol_in_rates[kinetics.id] > 0:
                                self.done_times[kinetics.id] = max_done_time_in + V/self.vol_in_rates[kinetics.id]
                            else:
                                self.done_times[kinetics.id] = max_done_time_in + V/self.volumetric_feed_rates[1]
                        #incoming conversion: Ca = Ca_0 * (1. - X) => X = 1. - Ca/Ca_0
                        conc_limiting, self.ready_flags[kinetics.id] = self.calc_conc(initial_conc = conc_in_sum,
                                                                                      initial_conversion = 1. - (conc_in_sum / self.conc0),
                                                                                      V=V,
                                                                                      reactor_type=kinetics.label,
                                                                                      k_eq=k_eq,
                                                                                      k=k,
                                                                                      id=kinetics.id)
                        self.conc_out_reactant[kinetics.id] = conc_limiting
                        self.conc_out_product[kinetics.id] = conc_product + conc_in_sum - conc_limiting   #taking into account the existing conc of products
                    else:#Do NOT output until ready
                        self.conc_out_reactant[kinetics.id] = 0.0
                        self.conc_out_product[kinetics.id] = 0.0
                        self.ready_flags[kinetics.id] = False
                    self.calc_outputs(kinetics.id, simulation_state, R)#now that this node has its outputs set, go on to the ones it outputs to

        if(not found):
            return


    async def calculate(self, simulation_state, graph):
        '''The actual simulation for number of objects specified by the protobuf '''
        #graph = graph # update the graph object when we get it (see controller.py)
        self.update_edge_list(graph)
        simulation_state = self.add_delete_protobuf_objects(simulation_state, graph)

        if (len(simulation_state.chemical_species) == 0):
            for i in range(len(self.chemical_species)):
                simulation_state.chemical_species.append(self.chemical_species[i])
        # for i in range(len(self.chemical_species)):
        #     simulation_state.chemical_species[i] = str(self.chemical_species[i])
        #TODO: Depending on how the selection of different reactions works within the code once updated within arc-board,
        #it may be necessary to add an "else" part of this if statement. This is because the simulation_state.chemical_species
        #will not be None, but it will potentially not be the correct chemical species either.

        if(not self.connected_to_source):
            self.start_plotting = False
            return simulation_state
        self.update_out_rates(0)#ONLY call this after update_edge_list() is done, and ONLY with id == 0

        if(self.reactor_number != len(simulation_state.kinetics)):#TODO: Find out why this is never(?) false...
            self.edge_list_changed = True

        if(self.edge_list_changed):  #reset simulation when edges change
            self.reactor_number = len(simulation_state.kinetics)
            self.start_time = self.graph_time
            self.edge_list_changed = False
            self.restart_plots = True
            self.start_plotting = True

        R = 8.314      # Universal gas constant (kJ/kmol K)


        self.time = self.graph_time - self.start_time
        self.conc_out_reactant, self.conc_out_product, self.conversion = {0:self.conc0}, {0:0}, {0:0}

        for kinetics in simulation_state.kinetics:
            self.conc_out_reactant[kinetics.id] = 0
            self.conc_out_product[kinetics.id] = 0
        self.calc_outputs(id = 0, simulation_state = simulation_state, R = R)
        #start with ID zero to look for nodes connected to source
        for kinetics in simulation_state.kinetics:
            flow_rate_limiting = self.conc_out_reactant[kinetics.id] * self.vol_in_rates[kinetics.id]
            flow_rate_out_product = self.conc_out_product[kinetics.id] * self.vol_in_rates[kinetics.id]  #taking into account the existing conc of products
            molar_flow = [flow_rate_limiting, flow_rate_out_product]
            if all([i == 0 for i in molar_flow]):
                mole_frac = [0*item for item in molar_flow]
            else:
                mole_frac = [item/sum(molar_flow) for item in molar_flow]
            for j in range(len(molar_flow)):
                kinetics.mole_fraction.append(float(mole_frac[j]))
                kinetics.molar_flow_rate.append(float(molar_flow[j]))
        simulation_state.time = self.time
        return simulation_state


    def cstr(self, initial_conc, initial_conversion, V, k_eq = 5, k = 0.1): #TODO: Need to create bool (or similar) to prevent CSTR from working when a gaseous reaction.
        '''Calculates concentrations for a first order, reversible reaction in a CSTR.
        Parameters
        ----------
        initial_conc : float
            Concentration of limiting reactant entering the reactor
        V : float
            Denotes the total volume of the reactor
        k_eq : float
            Denotes equilibrium concentration of the reaction
        k : float
            Denotes the reaction constant for the forward reaction
        Returns
        -------
        float
            Final concentration of the limiting reactant when it leaves the reactor
        '''

        cumulative_conversion = min(( initial_conversion + k * V) / (1. + (k + k/k_eq) * V), 1.) #since our volume flow rate is 1.0, we have tau = V. #TODO: generalize the (1-initial_conc), i.e. incoming conversion, for splits and any volumetric flow rate
        out_conc_lr = self.conc0*(1.0 - cumulative_conversion)
        ready = True #CSTR is instantaneous
        return (out_conc_lr, ready)

    def pfr(self, initial_conc, V, k_eq = 5, k = 0.1, done_time = None):
        '''Calculates concentrations for a first order, reversible reaction in a PFR.
        Note that incoming concentration only matters in finding final concentration due to mole balance and rate equation.
        Parameters
        ----------
        initial_conc : float
            Concentration of limiting reactant entering the reactor
        V : float
            Denotes the total volume of the reactor
        k_eq : float
            Denotes equilibrium concentration of the reaction
        k : float
            Denotes the reaction constant for the forward reaction
        Returns
        -------
        float
                Final concentration of the limiting reactant when it leaves the reactor
        '''

        if(done_time is None):
            done_time = V/self.volumetric_feed_rates[1]
        factor = 5.0 # 3.25
        time = min(done_time, self.time/factor) #divide by factor to arbitrarily accelerate display
        ready = False
        if(self.time/factor >= done_time):
            ready = True
        conversion = min(k_eq  / (k_eq + self.c/self.a) * (1. - math.exp( -time * k * (( self.c/self.a + k_eq ) / k_eq) ) ), 1.)
        out_conc_lr = initial_conc*(1.0 - conversion)
        return (out_conc_lr, ready)

    def pbr(self, initial_conc, R, k_eq = 5., k = 0.1): ##NOTE: This is actually for a BATCH REACTOR.
        '''Calculates concentrations for a first order, reversible reaction in a BR.
        Parameters
        ----------
        initial_conc : float
            Concentration of limiting reactant entering the reactor
        k : float
            Denotes the reaction constant for the forward reaction
        k_eq : float
            Denotes the equilibrium constant for the reaction
        Returns
        -------
        float
                Final concentration of the limiting reactant when it leaves the reactor
        '''
        t = self.time * 100 #making batch reactor instantaneous, seconds
        alpha = (1. + 1./k_eq) #for tidyness

        conversion =  1./alpha * ( 1. - math.exp(-alpha * k * t) )
        out_conc_lr = initial_conc * (1.0 - conversion)
        return (out_conc_lr, False) # Batch reactors never pump out

    # def pbr(self, initial_conc, V, k_eq = 5, k = 0.1):
    #     '''Calculates concentrations for a first order, reversible reaction in a PBR.
    #     Parameters
    #     ----------
    #     initial_conc : float
    #         Concentration of limiting reactant entering the reactor
    #     V : float
    #         Denotes the total volume of the reactor
    #     t: int
    #         Time in simulation at which concentration needs to be calculated
    #     k_eq : float
    #         Denotes equilibrium concentration of the reaction
    #     k : float
    #         Denotes the reaction constant for the forward reaction
    #     Returns
    #     -------
    #     float
    #             Final concentration of the limiting reactant when it leaves the reactor
    #     '''

    #     conversion = min((math.exp(V * k * (1 + self.c / self.a / k_eq)/self.vol_in_rates[kinetics.id]) + 1)/(1 + self.c / self.a / k_eq), 1)
    #     out_conc_lr = initial_conc * (1.0 - conversion)
    #     return (out_conc_lr, ready)

