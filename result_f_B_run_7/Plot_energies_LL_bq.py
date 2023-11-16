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
    a1.scatter(energies[:,0] * np.pi / E_c, energies[:,1], marker = "x")
    #a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    
    plt.show()
    
    return


def print_energies(parameters):
    """Print energies sector wise"""
    v_F = parameters[0,1]
    K = parameters[1,1]
    U_m = parameters[3,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1])
    
    print("-------------")
    print("K = " + str(K) + ", U_m = " + str(U_m))
    print("-------------")
    
    # plot
    fig = plt.figure()
    a1 =  fig.add_subplot(1,1,1)
    a1.set_title("Microscopic energies " , fontsize=30)
    a1.set_xlabel("k", fontsize=30)
    a1.set_ylabel("E", fontsize=30)
    #a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    
    for l in range(0, 2 * E_c + 1):
        total_momentum = (l - E_c) * 2 * np.pi / L
        
        # get energies 
        filename_energies_l = "energies/energies_J_1_" + str(l) + ".csv" 
        
        energies_l = np.genfromtxt(filename_energies_l, dtype=float, delimiter=' ')
        
        print("Current index l = " + str(l) + ", momentum p_l = " + str(total_momentum))
        print(energies_l)
        
        a1.scatter(total_momentum * np.ones(len(energies_l)), energies_l, marker = "x", color = "r")
    
    # calculate momenta
    k_vals = np.linspace(-E_c, E_c, 2 * E_c + 1) * 2 * np.pi / L
    
    # calculate expectation
    spectrum_expected = np.abs(k_vals) * (v_F / (2 * K))
    
    a1.plot(k_vals, spectrum_expected, color = "b", linestyle = "-",
            label = "Expected spectrum" ) 
    
    plt.show()
    
    


def main():      
    # Read in data
    
    parameters = np.genfromtxt("parameters.csv", delimiter=' ')
    
    
    print_energies(parameters)
    
    
    
    
 


if __name__ == '__main__':
    main()




	