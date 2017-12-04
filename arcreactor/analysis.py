import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time


class Analyzer:
    '''Controls plotting of different object properties'''
    def __init__(self):
        self.plot_number = 1
        self.xdata = []
        self.reactors = []
        self.reactor_number = 0
        self.start_time = 0

    @property
    def stream_names(self):
        return {'Reactor': ['plot']}

    def get_plot(self, name, simulation_state, start_plotting, restart_plots):
        if(restart_plots):
            self.xdata = []
            self.reactors = []
            self.reactor_number = 0
            self.start_time = simulation_state.time
        if(start_plotting):
            return self.plot_reactors(simulation_state)
        return None

    def plot_reactors(self,simulation_state):
        '''Plots reactor concentrations as received from simulation

        Parameters
        ----------
        simulation_state : protobuf object
                        Contains a time stamp and object properties

        Returns
        -------
        callable
                  Callable that retrieves contents of the output file

        '''
        if(len(simulation_state.kinetics) == 0):
            self.reactor_number = 0
            return None
        if(self.reactor_number != len(simulation_state.kinetics)):
            self.reactor_number = len(simulation_state.kinetics)
            self.start_time = simulation_state.time
            self.xdata = []
            self.reactors = []

        x = (simulation_state.time - self.start_time)*0.04      #display time in seconds(considering ~25fps)
        self.xdata.append(x)
        for i in range(len(simulation_state.kinetics)):
            y = simulation_state.kinetics[i].mole_fraction
            if(len(self.reactors) == i):#add one if we're out of space
                self.reactors.append([[] for j in range(len(y))])
            for j in range(len(y)):#reactors stores the concentrations of each species of reactant in each reactor over time
                self.reactors[i][j].append(y[j])#append this time step's molefrac to the respective reactor/species

        if(simulation_state.time % 25 == 0):
            fig, axes = plt.subplots(len(simulation_state.kinetics), 1,  sharex = True, sharey = True, squeeze=False) #overlay all the plots
            fig.text(0.5, 0.04, 'Time Elapsed (s)', ha='center', va='center')
            fig.text(0.01, 0.5, 'Concentration (mol/L)', ha='center', va='center', rotation='vertical')
            labels = ['A', 'B', 'C', 'D']
            colors = ['red', 'blue', 'green', 'purple']
            linestyles = ['-', ':', '-', ':']
            i = 0
            for i,ax in enumerate(axes[:,0]):
                #self.reactors[i] = [[] for j in range(len(y))]
                for j in range(len(y)):
                    #print(self.xdata, ydata[:,j])
                    ax.plot(self.xdata, self.reactors[i][j], color = colors[j], label = labels[j], ls=linestyles[j])

            plt.legend(loc='upper right')
            with io.BytesIO() as output:
                fig.savefig(output, format='jpg')
                plt.clf()
                return output.getvalue()