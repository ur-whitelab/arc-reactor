import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time


class Analyzer:

    def __init__(self, start_time):
        self.plot_number = 1
        self.start_time = start_time
        self.xdata = []
        self.r = []

    @property
    def stream_names(self):
        return {'Reactor': ['plot']}

    def get_plot(self, name, simulation_state):
        return self.plot_reactors(simulation_state)

    def plot_reactors(self,simulation_state):
        if(len(simulation_state.kinetics) == 0):
            return None

        fig, axes = plt.subplots(len(simulation_state.kinetics), sharex = True, sharey = True)
        #fig.axis([0,600, 0, 2])
        labels = ['C2H5COOCH3', 'H2O', 'CH3COOH', 'C2H5OH']
        colors = ['b', 'g', 'r', 'y']

        x = simulation_state.time - self.start_time          #time
        self.xdata.append(x)
        for i in range(len(simulation_state.kinetics)):
            y = simulation_state.kinetics[i].mole_fraction
            if(len(self.r) < (i+1)):
                self.r.append([])
            self.r[i].append(y)
            ydata = np.asarray(self.r[i])

            for j in range(len(y)):
                axes[i].plot(xdata,ydata[:,j], color = colors[j], label = labels[j])

        with io.BytesIO() as output:
            fig.savefig(output, format='jpg')
            plt.clf()
            return output.getvalue()
