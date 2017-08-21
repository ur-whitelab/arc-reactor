import io
import matplotlib.pyplot as plt
import numpy as np

class Analyzer:

    def __init__(self):
        self.plot_number = 1

    def get_plot(self, i, simulation_state):
        if i == 0:
            return self.plot_reactors(simulation_state)

    def plot_reactors(self,simulation_state):
        plt.figure()
        x = np.linspace(0,2*np.pi,100)
        y = np.sin(x + simulation_state.time / 20)
        plt.plot(x,y)
        output = io.BytesIO()
        plt.savefig(output, format='jpg')
        return output.getvalue()