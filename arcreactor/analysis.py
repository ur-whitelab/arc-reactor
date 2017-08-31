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

    def get_plot(self, i, simulation_state):
        if i == 0:
            return self.plot_reactors(simulation_state)

    def plot_reactors(self,simulation_state):
        fig, axes = plt.subplots(len(simulation_state.kinetics), sharex = True, sharey = True)
        #fig.axis([0,600, 0, 2])
        labels = ['C2H5COOCH3', 'H2O', 'CH3COOH', 'C2H5OH']
        colors = ['b', 'g', 'r', 'y']

        x = simulation_state.time - self.start_time          #time
        self.xdata.append(x)
        for i in range(len(simulation_state.kinetics)):
            y = simulation_state.kinetics[i].mole_fraction
            if(not self.r[i]):
                self.r.append([])
            self.r[i].append(y)
            ydata = np.asarray(self.r[i])

            for j in range(len(y)):
                axes[i].plot(xdata,ydata[:,j], color = colors[j], label = labels[j])

        axes[2].set_xlabel('Time (seconds)')
        axes[1].set_ylabel('Concentration (moles/dm3)')
        fig.legend()
        output = io.BytesIO()
        fig.savefig(output, format='jpg')
        return output.getvalue()