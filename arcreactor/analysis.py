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
        self.r = []
        self.reactor_number = 0
        self.start_time = 0

    @property
    def stream_names(self):
        return {'Reactor': ['plot']}

    def get_plot(self, name, simulation_state):
        return self.plot_reactors(simulation_state)

    def plot_reactors(self,simulation_state):
        '''Plots reactor concentrations as received from simulation

        Parameters
        ----------
        simulation_state : protobuf object
                        Protobuf object contains a time stamp and object properties

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
            self.r = []
        fig, axes = plt.subplots(len(simulation_state.kinetics), 1,  sharex = True, sharey = True, squeeze=False)
        labels = ['C2H5COOCH3', 'H2O', 'CH3COOH', 'C2H5OH']
        colors = ['b', 'g', 'r', 'y']
        x = (simulation_state.time - self.start_time)*0.04      #display time in seconds(considering ~25fps)
        self.xdata.append(x)
        i = 0
        for i,ax in enumerate(axes[:,0]):
            y = simulation_state.kinetics[i].mole_fraction
            if(len(self.r) == i):
                self.r.append([])
            self.r[i] = y
            ydata = np.array(self.r[i])
            for j in range(len(y)):
                #print(self.xdata, ydata[:,j])
                ax.plot(self.xdata, ydata[:,j], color = colors[j], label = labels[j])
            plt.legend()

        with io.BytesIO() as output:
            fig.savefig(output, format='jpg')
            plt.clf()
            return output.getvalue()
