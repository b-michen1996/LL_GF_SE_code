import csv
import numpy as np
import matplotlib.animation as animation
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib
import kwant



def plot_energies(energies, parameters):
    """Plot syndromes data"""
    # get some parameters
    L = parameters[5,1]
    E_c = int(parameters[7,1])
    p_c = int(parameters[8,1])    
    omega = int(parameters[9,1])    
    
    # calculate momenta
    k_vals = np.linspace(-E_c, E_c, 2 * E_c + 1) * np.pi / E_c
        
    # plot
    fig = plt.figure()
    a1 =  fig.add_subplot(1,1,1)
    a1.set_title("Microscopic energies " , fontsize=30)
    a1.set_xlabel("k", fontsize=30)
    a1.set_ylabel("E", fontsize=30)
    a1.scatter(energies[:,0], energies[:,1], marker = "x")
    a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    
    plt.show()
    
    return
    


def main():      
    # Read in data
    energies = np.genfromtxt("energies.csv", dtype=complex, delimiter=' ')
    parameters = np.genfromtxt("parameters.csv", delimiter=' ')
    
    plot_energies(energies, parameters)
    
    
    
    
    
    
 


if __name__ == '__main__':
    main()




	