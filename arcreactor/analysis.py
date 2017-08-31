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
        self.r1 = []
        self.r2 = []
        self.r3 = []

    def get_plot(self, i, simulation_state):
        if i == 0:
            return self.plot_reactors(simulation_state)

    def plot_reactors(self,simulation_state):
        fig, axes = plt.subplots(3, sharex = True, sharey = True)
        fig.axis([0,600, 0, 2])
        labels = ['C2H5COOCH3', 'H2O', 'CH3COOH', 'C2H5OH']
        colors = ['b', 'g', 'r', 'y']
        
        x = self.start_time - simulation_state.time          #time
        self.xdata.append(x)
        for i in range(3):
            y = simulation_state.kinetics[i].mole_fraction
            if i == 0:
                self.r1.append(y)
                ydata = np.asarray(self.r1)
            if i ==1:
                self.r2.append(y)
                ydata = np.asarray(self.r2)
            if i == 2:
                self.r3.append(y)
                ydata = np.asarray(self.r3)
            
            for j in range(4):
                axes[i].plot(xdata,ydata[:,j], color = colors[j], label = labels[j])

        axes[2].set_xlabel('Time (seconds)')
        axes[1].set_ylabel('Concentration (moles/dm3)')
        fig.legend()
        output = io.BytesIO()
        fig.savefig(output, format='jpg')
        return output.getvalue()