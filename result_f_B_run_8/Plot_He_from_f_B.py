import csv
import numpy as np
import matplotlib.animation as animation
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib
import kwant



def plot_spectrum(parameters, omega, beta, eta):
    """Plot syndromes data"""
    # get some parameters
    K = parameters[1,1]
    alpha = parameters[2,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
    
    # calculate GF
    GF = calc_GF_from_f_B(parameters, omega, beta, eta)
    
    # calculate spectrum    
    spectrum = calc_spectrum_from_GF(GF, omega)
    
    # calculate momenta
    k_vals = np.linspace(-p_c, p_c, 2 * p_c + 1) * np.pi / p_c
    
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
    a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    a1.legend(loc='upper center', prop={'size': 12})
    # change color of legend to black
    leg = a1.get_legend()
    leg.legendHandles[0].set_color('black')
    leg.legendHandles[1].set_color('black') 
    
    plt.show()
    
    return
    

def calc_spectrum_from_GF(GF, omega):
    """ Calculate spectrum of H_e from GF for given omega."""    
    dimensions = GF.shape
    
    spectrum= 1j * np.zeros((dimensions[0], 2))
    
    # loop over all momenta
    for l in range(0, dimensions[0]):
        # get GF for current momentum
        GF_k = 1j * np.zeros((2,2))
        GF_k[0,0] = GF[l, 0]
        GF_k[0,1] = GF[l, 1]
        GF_k[1,0] = GF[l, 2]
        GF_k[1,1] = GF[l, 3]
        
        # get effective Hamiltonian 
        H_e_k = 1j * np.zeros((2,2))
        
        
        try:
            H_e_k = omega * np.zeros((2,2)) - np.linalg.inv(GF_k)
        except:
            print("Warning: singular GF at l = " + str(l) + "!")
            print(GF_k)
        
        # calculate and store eigenvalues
        eigenvalues, eigenvectors = np.linalg.eig(H_e_k)
        
        spectrum[l, 0] = eigenvalues[0]
        spectrum[l, 1] = eigenvalues[1]
        
    return spectrum

            
def calc_GF_from_f_B(parameters, omega, beta, eta):
    """ Read in matrix-representation of bosonic factors form file and calculate
        GF for given omega, beta, eta."""
    # get some parameters
    K = parameters[1,1]
    alpha = parameters[2,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
    
    GF = 1j * np.zeros((int(2 * p_c + 1), 4))
    
    # calculate partition sum and global prefactor
    Z = 0
    for l in range(0, 2 * E_c + 1): 
        # get energies of sector
        filename_energies_l = "energies/energies_J_1_" + str(l) + ".csv" 
        energies_l = np.genfromtxt(filename_energies_l, dtype=float, delimiter=' ')
        
        # add to sum 
        Z += np.sum(np.exp(- beta * energies_l))
    prefactor = (2 * np.pi * alpha / L) ** ((1-K)**2 / (2*K)) / Z

    # loop over all momenta from -p_c to p_c
    for l in range(0, 2 * p_c + 1):
        k = l - p_c
        
        # loop over all momentum sectors with non-vanishing matrix elements of f_B        
        for l2 in range(max(0, k) , 2 * E_c + 1):
            l1 = l2 - k               
            if l1 >  2 * E_c:
                break
            
            # get matrix representation of f_B in eigenbasis
            filename_f_B_R_l1_l2 = "f_B_R_" + str(l) + "/f_B_R_" + str(l1) + "_" + str(l2) + ".csv"
            filename_f_B_L_l1_l2 =  "f_B_L_" + str(l) + "/f_B_L_" + str(l1) + "_" + str(l2) + ".csv"
            
            f_B_R = np.genfromtxt(filename_f_B_R_l1_l2, dtype=complex, delimiter=' ')
            f_B_L = np.genfromtxt(filename_f_B_L_l1_l2, dtype=complex, delimiter=' ')
            
            # get energies 
            filename_energies_l1 = "energies/energies_J_1_" + str(l1) + ".csv"
            filename_energies_l2 = "energies/energies_J_1_" + str(l2) + ".csv"
            
            energies_l1 = np.genfromtxt(filename_energies_l1, dtype=float, delimiter=' ')
            energies_l2 = np.genfromtxt(filename_energies_l2, dtype=float, delimiter=' ')

            # carry out summation for sector
            G_RR = summation_sector(f_B_R, f_B_R, energies_l1, energies_l2, omega, beta, eta)
            G_RL = summation_sector(f_B_R, f_B_L, energies_l1, energies_l2, omega, beta, eta)
            G_LR = summation_sector(f_B_L, f_B_R, energies_l1, energies_l2, omega, beta, eta)
            G_LL = summation_sector(f_B_L, f_B_L, energies_l1, energies_l2, omega, beta, eta)
            
            if G_LL == 0:
                print("Summation is zero for l1 = ", l1, ", l2 = ", l2)
                
                
            # print(summation_sector(f_B_R, f_B_R, energies_l1, energies_l2, omega, beta, eta))
            # write to GF array
            GF[l, :] = prefactor * np.array([G_RR, G_RL, G_LR, G_LL])
    return GF
                    

def summation_sector(f_B_r1, f_B_r2, energies_l1, energies_l2, omega, beta, eta):
    """"Carry out the summation from the Kaellen-Lehmann representation for 
    the given parameters"""
    dim_l1 = len(energies_l1)
    dim_l2 = len(energies_l2)
    
    result = 0j
        
    for j1 in range (0, dim_l1):
        for j2 in range (0, dim_l2):
            result = result + f_B_r1[j1, j2] * np.conj(f_B_r2[j1, j2]) / (omega 
                            + energies_l1[j1] - energies_l2[j2] + 1j * eta)

    return result
    

def main():      
    # Read in data
    parameters = np.genfromtxt("parameters.csv", delimiter=' ')
    
    omega= 0
    beta = 1000
    eta = 0.1
    
    plot_spectrum(parameters, omega, beta, eta)


if __name__ == '__main__':
    main()




	