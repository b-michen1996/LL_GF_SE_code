import csv
import numpy as np
import matplotlib.animation as animation
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib
import kwant



def plot_spectrum(G_RR, G_RL, G_LR, G_LL, parameters):
    """Plot syndromes data"""
    # get some parameters
    L = parameters[5,1]
    E_c = int(parameters[7,1])
    p_c = int(parameters[8,1])    
    omega = parameters[9,1]

    
    # calculate spectrum    
    spectrum = calc_spectrum(G_RR, G_RL, G_LR, G_LL, p_c, omega)
    
    # calculate momenta
    k_vals = np.linspace(-p_c, p_c, 2 * p_c + 1) * 2 * np.pi / L
    
    # sort energies 
    energies_re = np.sort(spectrum.real, 1)
    energies_im = np.sort(spectrum.imag, 1)
    
    # plot
    fig = plt.figure()
    a1 =  fig.add_subplot(1,1,1)
    a1.set_title("Spectrum of $H_{eff}$ " , fontsize=30)
    a1.set_xlabel("k", fontsize=30)
    a1.set_ylabel("E", fontsize=30)
    a1.plot(k_vals, energies_re[:,0], 'r-', label = "Re(E)") 
    a1.plot(k_vals, energies_im[:,0],'r--', label = "Im(E)")
    a1.plot(k_vals, energies_re[:,1], 'b-')
    a1.plot(k_vals, energies_im[:,1],'b--')  
    #a1.set_xlim([-np.pi/2, np.pi/2])
    #a1.set_ylim([-0.1, 0.1])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    a1.legend(loc='upper center', prop={'size': 12})
    # change color of legend to black
    leg = a1.get_legend()
    leg.legendHandles[0].set_color('black')
    leg.legendHandles[1].set_color('black') 
    
    plt.show()
    
    return
    
def calc_spectrum(G_RR, G_RL, G_LR, G_LL, p_c, omega):
    """Calculate complex spectrum of effective Hamiltonian from GF"""
    spectrum= 1j * np.zeros((int(2 * p_c + 1), 2))
    
    for l in range(0, 2 * p_c + 1):
        GF_k = 1j * np.zeros((2,2))
        GF_k[0,0] = G_RR[l]
        GF_k[0,1] = G_RL[l]
        GF_k[1,0] = G_LR[l]
        GF_k[1,1] = G_LL[l]
        
        H_e_k = omega * np.identity(2)
        
        try:
            H_e_k = H_e_k - np.linalg.inv(GF_k)
        except:
            print("Warning: singular GF at k = " + str(l - p_c) + "!")
                 
        eigenvalues, eigenvectors = np.linalg.eig(H_e_k)
        
        spectrum[l, 0] = eigenvalues[0]
        spectrum[l, 1] = eigenvalues[1]
        
        print("l = ", l, ", k = ", l - p_c, "E = ", spectrum[l, :])
    
    return spectrum


def main():      
    # Read in data
    G_RR = np.genfromtxt("G_RR.csv", dtype=complex, delimiter=' ')
    G_RL = np.genfromtxt("G_RL.csv", dtype=complex, delimiter=' ')
    G_LR = np.genfromtxt("G_LR.csv", dtype=complex, delimiter=' ')
    G_LL = np.genfromtxt("G_LL.csv", dtype=complex, delimiter=' ')

    parameters = np.genfromtxt("parameters.csv", delimiter=' ')
    
    plot_spectrum(G_RR, G_RL, G_LR, G_LL, parameters)
    
    
    
    
    
    
 


if __name__ == '__main__':
    main()




	