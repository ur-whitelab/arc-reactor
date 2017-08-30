import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time

plt.ion()

class Analyzer:

    def __init__(self, start_time):
        self.plot_number = 1
        self.start_time = start_time

    def get_plot(self, i, simulation_state):
        if i == 0:
            return self.plot_reactors(simulation_state)

    def plot_reactors(self,simulation_state):
        x = self.start_time - simulation_state.time          #time
        #y = simulation_state.kinetics[0].mole_fraction
        labels = ['C2H5COOCH3', 'H2O', 'CH3COOH', 'C2H5OH']
        colors = ['b', 'g', 'r', 'y']
        for i in range(3):
            y = simulation_state.kinetics[i].mole_fraction
            axes[i].set_title('Reactor {}'.format(i+1))
            for j in range(3):
                axes[i].plot(x,y[j], color = colors[i], label = labels[j])
        axes[2].set_xlabel('Time (seconds)')
        axes[1].set_ylabel('Concentration (moles/dm3)')
        fig.legend()
        output = io.BytesIO()
        fig.savefig(output, format='jpg')
        return output.getvalue()

