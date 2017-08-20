import io
import matplotlib.pyplot as plt
import numpy as np

def plot_reactors(simulation_state):
    plt.figure()
    x = np.linspace(0,2*np.pi,100)
    y = np.sin(x + simulation_state.time / 20)
    plt.plot(x,y)
    output = io.BytesIO()
    plt.savefig(output, format='jpg')
    return output.getvalue()