import numpy as np
import asyncio
import datetime as dt
import time
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Hydrolysis of ester - A psuedo first order chemical reaction
# All reactors are equal size and operate at the same temperature
#All numbers are arbitrary

no_of_reactors = 3
reactor_volumes = [2,2,2]		#dm3
reactants = 2	
A = 'CH3COOC2H5'
B = 'H2O'
C = 'CH3COOH'
D = 'C2H5OH'
k = 0.01
feed_volumetric_rates = [0.5,0.5]		#dm3/s
feed_rate = [1,1]					#mol/s

def rxn_d(c,t,k=0.01):
	return -k*c

def calculate(start_time, reactor_volumes, feed_volumetric_rates, feed_rate):
	'''Calculates all the kinetics of the reactors'''
	c0 = np.asarray(feed_rate)/np.asarray(feed_volumetric_rates)		#mol/dm3
	tau = reactor_volumes[0]/feed_volumetric_rates[0]
	t_int = np.linspace(0, 3600, 3600000)								#simulate for one hour 
	ca1 = odeint(rxn_d, c0[0], t_int) 
	x1 = 1 - ca1/c0[0]
	x = int(dt.datetime.now().microsecond/1000 - start_time.microsecond/1000)

	# The following steps for calculating reactor concentrations can be looped later
	c1 = np.asarray([ca1[x], ca1[x], c0[0] - ca1[x], c0[0] - ca1[x]])
	mf1 = c1/np.sum(c1)
	ca2 = ca1/(1 + tau*k)
	c2 = np.asarray([ca2[x], ca2[x], c0[0] - ca2[x], c0[0] - ca2[x]])
	mf2 = c2/np.sum(c2)
	ca3 = ca2/(1 + tau*k)
	c3 = np.asarray([ca3[x], ca3[x], c0[0] - ca3[x], c0[0] - ca3[x]])
	mf3 = c3/np.sum(c3)
	return mf1, mf2, mf3





