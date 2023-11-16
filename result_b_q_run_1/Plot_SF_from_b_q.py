import csv
import numpy as np
import matplotlib.animation as animation
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib
import kwant
import time

# enable latex
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "times"
})

    
def plot_SF_bq(parameters, n_omega, beta, eta):
    """Plot spectrum of self-energy obtained from b_q correlator"""
    # get some parameters
    v_F = parameters[0,1]
    K = parameters[1,1]
    alpha = parameters[2,1]
    U_m = parameters[3,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
    
    title = (r"A(k, $\omega$) for $\beta$ = " 
             + str(beta) + r", $\eta$ = " + str(eta)) 
    
    # plot
    fig = plt.figure(figsize=(24, 12), layout = "tight")
    a1 =  fig.add_subplot(1,1,1)
    a1.set_title(title , fontsize=30)
    a1.set_xlabel("k", fontsize=30)
    a1.set_ylabel("E", fontsize=30)
    
    # calculate momenta
    k_vals = np.linspace(-p_c, p_c, 2 * p_c + 1) * 2 * np.pi / L
    
    # calculate expectation
    spectrum_expected = np.abs(k_vals) * (v_F / (2 * K))
    
    omega_vals =  (1.2 / n_omega) * max(spectrum_expected) * np.linspace(-1, n_omega, n_omega) 
        
    SF = np.zeros((n_omega, 2 * p_c))
    
    K =  np.zeros((n_omega, 2 * p_c))
    W = np.zeros((n_omega, 2 * p_c))
    
    for j in range(n_omega):
        omega = omega_vals[j]
        # calculate GF
        GF = calc_GF_from_bq(parameters, omega, beta, eta)
        SF[j,:] = -(1 / np.pi) * GF.imag
        W[j:] = omega * np.ones(2 * p_c)
    
    for j in range(0 , p_c):
        K[:,j] = (j - p_c) * 2 * np.pi / L 
    
    for j in range(p_c , 2 * p_c):
        K[:,j] = (j - p_c + 1) * 2 * np.pi / L 

    cax = a1.contourf(K, W, SF, levels = 100, cmap=plt.cm.magma_r)
    
    a1.plot(k_vals, spectrum_expected, color = "g", linestyle = "-",
            label = "Expected spectrum" ) 
    
    """
    a1.plot(k_vals, spectrum.real, color = "r", linestyle = "-",
            label = r"$\omega$ = " + str(omega) + r", $\beta$ = " + str(beta) + r", $\eta$ = " + str(eta) ) 
    a1.plot(k_vals, spectrum.imag, color = "r",  linestyle = "--")
    
    a1.plot(k_vals, spectrum_expected, color = "b", linestyle = "-",
            label = "Expected spectrum" ) 
    """
    #a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    # change color of legend to black
    leg = a1.get_legend()
    fig.colorbar(cax)
    fig.savefig("SF_beta_" + str(beta) + "_eta_" + str(eta) + ".png",  bbox_inches='tight')
    
    plt.show()
    
    return


def calc_spectrum_from_GF(GF, omega, p_c):
    """ Calculate spectrum of H_e from GF for given omega."""        
    spectrum= 1j * np.zeros(2 * p_c + 1)
    
    spectrum[:p_c] = omega * np.ones(p_c) - 1 / GF[:p_c]
    spectrum[p_c] = omega 
    spectrum[p_c + 1:] = omega * np.ones(p_c) - 1 / GF[p_c:]
        
    return spectrum


def calc_GF_from_bq(parameters, omega, beta, eta):
    """ Read in matrix-representation of b_q from file and calculate
        GF for given omega, beta, eta."""
    # get some parameters
    K = parameters[1,1]
    alpha = parameters[2,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
    
    GF = 1j * np.zeros(int(2 * p_c))
    
    # calculate partition sum and global prefactor
    Z = 0
    for l in range(0, 2 * E_c + 1): 
        # get energies of sector
        filename_energies_l = "energies/energies_J_1_" + str(l) + ".csv" 
        energies_l = np.genfromtxt(filename_energies_l, dtype=float, delimiter=' ')
                
        # add to sum 
        Z += np.sum(np.exp(- beta * energies_l))
        
    print("partition sum Z = ", Z)
    prefactor = 1 / Z

    # loop over all momenta from -p_c to p_c
    for l in range(0, 2 * p_c):
        q = l - p_c + 1
        if (l < p_c):
            q = l - p_c
        
        # loop over all momentum sectors with non-vanishing matrix elements of f_B        
        for l2 in range(max(0, q) , 2 * E_c + 1):
            l1 = l2 - q               
            if l1 >  2 * E_c:
                break
            tic =time.perf_counter()
            # get matrix representation of f_B in eigenbasis
            filename_b_ql_l1_l2 = "b_" + str(l) + "/b_" + str(l1) + "_" + str(l2) + ".csv"
            
            b_ql_l1_l2 = np.genfromtxt(filename_b_ql_l1_l2, dtype=complex, delimiter=' ')
            
            # get energies 
            filename_energies_l1 = "energies/energies_J_1_" + str(l1) + ".csv"
            filename_energies_l2 = "energies/energies_J_1_" + str(l2) + ".csv"
            
            energies_l1_data = np.genfromtxt(filename_energies_l1, dtype=float, delimiter=' ')
            energies_l2_data = np.genfromtxt(filename_energies_l2, dtype=float, delimiter=' ')
            

                        
            d1 = len(energies_l1_data)
            d2 = len(energies_l2_data)
            # create matrix
            energies_l1 = np.ones((d1, d2))
            energies_l2 = np.ones((d1, d2))
            
            for j in range(0, d1):
                energies_l1[j, :] = energies_l1_data[j] * np.ones(d2)
                
            for j in range(0, d2):
                energies_l2[:, j] = energies_l2_data[j] * np.ones(d1)
            
            toc =time.perf_counter()
            #print("time for reading in data: ", toc-tic)
            
            tic =time.perf_counter()
            # carry out summation for sector            
            G_contribution = summation_sector_bq(b_ql_l1_l2, energies_l1, energies_l2, omega, beta, eta)

            toc =time.perf_counter()
            #print("time for summation: ", toc-tic)
            
            # print(summation_sector(f_B_R, f_B_R, energies_l1, energies_l2, omega, beta, eta))
            # write to GF array
            GF[l] += prefactor * G_contribution
    return GF


def summation_sector_bq(b_ql_l1_l2, energies_l1, energies_l2, omega, beta, eta):
    """"Carry out the summation from the Kaellen-Lehmann representation for 
    the given parameters"""
    
    factor_matrix = np.divide(np.exp(-beta * energies_l1) - np.exp(-beta * energies_l2),
    (omega + 1j *eta) * np.ones(energies_l1.shape) + energies_l1  - energies_l2)
        
    result = np.sum(b_ql_l1_l2 * np.conj(b_ql_l1_l2) * factor_matrix)
    
    return result
    

def main():      
    # Read in data
    parameters = np.genfromtxt("parameters.csv", delimiter=' ')
    
    beta = 100
    eta = 0.001
    
    n_omega = 100

    plot_SF_bq(parameters, n_omega, beta, eta)


if __name__ == '__main__':
    main()




	