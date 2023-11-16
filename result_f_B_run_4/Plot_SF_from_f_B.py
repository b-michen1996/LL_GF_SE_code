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

    
def plot_SF_f_B_R(parameters, n_omega, beta, eta):
    """Plot spectrum of self-energy obtained from f_B_R correlator"""
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
    spectrum_expected_R = k_vals * (v_F / (2 * K))
    spectrum_expected_L = -k_vals * (v_F / (2 * K))
    
    omega_vals =  1.2 * np.linspace(min(spectrum_expected_R), 
                                    max(spectrum_expected_R) , n_omega) 
        
    SF = np.zeros((n_omega, 2 * p_c + 1))
    
    K =  np.zeros((n_omega, 2 * p_c + 1))
    W = np.zeros((n_omega, 2 * p_c + 1))
    
    for j in range(n_omega):
        omega = omega_vals[j]
        # calculate GF
        GF = calc_GF_from_f_B(parameters, omega, beta, eta)
        SF[j,:] = -(1 / np.pi) * (GF[:, 0].imag + GF[:, 3].imag)
        W[j,:] = omega * np.ones(2 * p_c + 1)
        K[j,:] = k_vals
    


    cax = a1.contourf(K, W, SF, levels = 100, cmap=plt.cm.magma_r)
    
    a1.plot(k_vals, spectrum_expected_R, color = "g", linestyle = "-",
            label = "Expected spectrum" ) 
    a1.plot(k_vals, spectrum_expected_L, color = "g", linestyle = "-") 
    
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
        
    print("partition sum Z = ", Z)
    prefactor = (2 * np.pi * alpha / L) ** ((1-K)**2 / (2*K)) / Z
    
    check_f_B_R = 0
    check_f_B_L = 0

    # loop over all momenta from -p_c to p_c
    for l in range(0, 2 * p_c + 1):
        k = l - p_c
        
        # loop over all momentum sectors with non-vanishing matrix elements of f_B        
        for l2 in range(max(0, k) , 2 * E_c + 1):
            l1 = l2 - k               
            if l1 >  2 * E_c:
                break
            tic =time.perf_counter()
            # get matrix representation of f_B in eigenbasis
            filename_f_B_R_l1_l2 = "f_B_R_" + str(l) + "/f_B_R_" + str(l1) + "_" + str(l2) + ".csv"
            filename_f_B_L_l1_l2 =  "f_B_L_" + str(l) + "/f_B_L_" + str(l1) + "_" + str(l2) + ".csv"
            
            f_B_R = np.genfromtxt(filename_f_B_R_l1_l2, dtype=complex, delimiter=' ')
            f_B_L = np.genfromtxt(filename_f_B_L_l1_l2, dtype=complex, delimiter=' ')
            
            check_f_B_R += np.linalg.norm(f_B_R)
            check_f_B_L += np.linalg.norm(f_B_L)
            
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
            G_RR = summation_sector(f_B_R, f_B_R, energies_l1, energies_l2, omega, beta, eta)
            G_RL = 1j * summation_sector(f_B_R, f_B_L, energies_l1, energies_l2, omega, beta, eta)
            G_LR = -1j * summation_sector(f_B_L, f_B_R, energies_l1, energies_l2, omega, beta, eta)
            G_LL = summation_sector(f_B_L, f_B_L, energies_l1, energies_l2, omega, beta, eta)

            toc =time.perf_counter()
            #print("time for summation: ", toc-tic)
            
            # print(summation_sector(f_B_R, f_B_R, energies_l1, energies_l2, omega, beta, eta))
            # write to GF array
            GF[l, :] += prefactor * np.array([G_RR, G_RL, G_LR, G_LL])
    
    return GF


def summation_sector(f_B_r1, f_B_r2, energies_l1, energies_l2, omega, beta, eta):
    """"Carry out the summation from the Kaellen-Lehmann representation for 
    the given parameters"""
    
    factor_matrix = np.divide(np.exp(-beta * energies_l1) + np.exp(-beta * energies_l2),
    (omega + 1j *eta) * np.ones(energies_l1.shape) + energies_l1  - energies_l2)
        
    result = np.sum(f_B_r1 * np.conj(f_B_r2) * factor_matrix)
    
    return result


def main():      
    # Read in data
    parameters = np.genfromtxt("parameters.csv", delimiter=' ')
    
    beta = 150
    eta = 0.001
    
    n_omega = 200

    plot_SF_f_B_R(parameters, n_omega, beta, eta)


if __name__ == '__main__':
    main()




	