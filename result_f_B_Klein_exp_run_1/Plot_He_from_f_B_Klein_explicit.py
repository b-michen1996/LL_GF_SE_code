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

def plot_spectrum(parameters, omega, beta, eta):
    """Plot syndromes data"""
    # get some parameters
    v_F = parameters[0,1]
    K = parameters[1,1]
    alpha = parameters[2,1]
    U_m = parameters[3,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
    
    beta = np.array(beta)
    eta = np.array(eta)    
    
    title = (r"$v_F$ = " + str(v_F) +  r", K = " + str(K) +  r", $U_m$ = " + str(U_m) + r", $\gamma$ = " + str(alpha)
    + r", L = " + str(L)) + r", $E_c$ = " + str(E_c) + r", $p_c$ =" + str(p_c) 
    # plot
    fig = plt.figure()
    a1 =  fig.add_subplot(1,1,1)
    a1.set_title(title , fontsize=30)
    a1.set_xlabel("k", fontsize=30)
    a1.set_ylabel("E", fontsize=30)
    
    for b in beta:
        for e in eta:
            # calculate GF
            GF = calc_GF_from_f_B(parameters, omega, b, e)
    
            # calculate spectrum    
            spectrum = calc_spectrum_from_GF(GF, omega)
    
            # calculate momenta
            k_vals = np.linspace(-p_c, p_c, 2 * p_c + 1) * 2 * np.pi / L
            # sort energies 
            energies_re = np.sort(spectrum.real, 1)
            energies_im = np.sort(spectrum.imag, 1) 
            
            color_curr = next(a1._get_lines.prop_cycler)['color']
            label_string =  r", $\beta$ = " + str(b) + r", $\eta$ =" + str(e) 
            
            a1.plot(k_vals, energies_re[:,0], color = color_curr, linestyle = "-",
                    label = label_string ) 
            a1.plot(k_vals, energies_im[:,0], color = color_curr,  linestyle = "--")
            a1.plot(k_vals, energies_re[:,1], color = color_curr,  linestyle = "-")
            a1.plot(k_vals, energies_im[:,1], color = color_curr,  linestyle = "--")  
    
    # calculate expectation
    spectrum_expected_R = k_vals * (v_F / (2 * K))
    spectrum_expected_L = -k_vals * (v_F / (2 * K))
    
    a1.plot(k_vals, spectrum_expected_R, color = "r", linestyle = "dashdot",
            label = "Expected spectrum" ) 
    a1.plot(k_vals, spectrum_expected_L, color = "r", linestyle = "dashdot") 

    #a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    a1.legend(loc='upper center', prop={'size': 12})
    # change color of legend to black
    leg = a1.get_legend()
    
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
            # print(GF_k)
        
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
        filename_energies_l_J_1 = "energies/energies_J_1_" + str(l) + ".csv" 
        energies_l_J_1 = np.genfromtxt(filename_energies_l_J_1, dtype=float, delimiter=' ')
        
        filename_energies_l_J_m1 = "energies/energies_J_m1_" + str(l) + ".csv" 
        energies_l_J_m1 = np.genfromtxt(filename_energies_l_J_m1, dtype=float, delimiter=' ')
                
        # add to sum 
        Z += np.sum(np.exp(- beta * energies_l_J_1)) + np.sum(np.exp(- beta * energies_l_J_m1))
        
    print("partition sum Z = ", Z)
    prefactor = (2 * np.pi * alpha / L) ** ((1-K)**2 / (2*K)) / Z

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
            filename_f_B_R_l1_l2_t1 = "f_B_R_" + str(l) + "/f_B_R_t1_" + str(l1) + "_" + str(l2) + ".csv"
            filename_f_B_L_l1_l2_t1 =  "f_B_L_" + str(l) + "/f_B_L_t1_" + str(l1) + "_" + str(l2) + ".csv"
            filename_f_B_R_l1_l2_t2 = "f_B_R_" + str(l) + "/f_B_R_t2_" + str(l1) + "_" + str(l2) + ".csv"
            filename_f_B_L_l1_l2_t2 =  "f_B_L_" + str(l) + "/f_B_L_t2_" + str(l1) + "_" + str(l2) + ".csv"
            
            f_B_R_t1 = np.genfromtxt(filename_f_B_R_l1_l2_t1, dtype=complex, delimiter=' ')
            f_B_L_t1 = np.genfromtxt(filename_f_B_L_l1_l2_t1, dtype=complex, delimiter=' ')
            f_B_R_t2 = np.genfromtxt(filename_f_B_R_l1_l2_t2, dtype=complex, delimiter=' ')
            f_B_L_t2 = np.genfromtxt(filename_f_B_L_l1_l2_t2, dtype=complex, delimiter=' ')
            
            # get energies 
            filename_energies_l1_t1 = "energies/energies_J_1_" + str(l1) + ".csv"
            filename_energies_l2_t1 = "energies/energies_J_m1_" + str(l2) + ".csv"
            filename_energies_l1_t2 = "energies/energies_J_m1_" + str(l1) + ".csv"
            filename_energies_l2_t2 = "energies/energies_J_1_" + str(l2) + ".csv"
            
            energies_l1_data_t1 = np.genfromtxt(filename_energies_l1_t1, dtype=float, delimiter=' ')
            energies_l2_data_t1 = np.genfromtxt(filename_energies_l2_t1, dtype=float, delimiter=' ')
            energies_l1_data_t2 = np.genfromtxt(filename_energies_l1_t2, dtype=float, delimiter=' ')
            energies_l2_data_t2 = np.genfromtxt(filename_energies_l2_t2, dtype=float, delimiter=' ')
            
            d1 = len(energies_l1_data_t1)
            d2 = len(energies_l2_data_t1)
            
            # create matrix
            energies_l1_t1 = np.ones((d1, d2))
            energies_l2_t1 = np.ones((d1, d2))
            energies_l1_t2 = np.ones((d1, d2))
            energies_l2_t2 = np.ones((d1, d2))
            
            for j in range(0, d1):
                energies_l1_t1[j, :] = energies_l1_data_t1[j] * np.ones(d2)
                energies_l1_t2[j, :] = energies_l1_data_t2[j] * np.ones(d2)
                
            for j in range(0, d2):
                energies_l2_t1[:, j] = energies_l2_data_t1[j] * np.ones(d1)
                energies_l2_t2[:, j] = energies_l2_data_t2[j] * np.ones(d1)
            
            toc =time.perf_counter()
            #print("time for reading in data: ", toc-tic)
            
            tic =time.perf_counter()
            # carry out summation for sector            
            G_RR_t1 = summation_sector(f_B_R_t1, f_B_R_t1, energies_l1_t1, energies_l2_t1, omega, beta, eta)
            G_RL_t1 = 1j * summation_sector(f_B_R_t1, f_B_L_t1, energies_l1_t1, energies_l2_t1, omega, beta, eta)
            G_LR_t1 = -1j * summation_sector(f_B_L_t1, f_B_R_t1, energies_l1_t1, energies_l2_t1, omega, beta, eta)
            G_LL_t1 = summation_sector(f_B_L_t1, f_B_L_t1, energies_l1_t1, energies_l2_t1, omega, beta, eta)
            
            G_RR_t2 = summation_sector(f_B_R_t2, f_B_R_t2, energies_l1_t2, energies_l2_t2, omega, beta, eta)
            G_RL_t2 = -1j * summation_sector(f_B_R_t2, f_B_L_t2, energies_l1_t2, energies_l2_t2, omega, beta, eta)
            G_LR_t2 = 1j * summation_sector(f_B_L_t2, f_B_R_t2, energies_l1_t2, energies_l2_t2, omega, beta, eta)
            G_LL_t2 = summation_sector(f_B_L_t2, f_B_L_t2, energies_l1_t2, energies_l2_t2, omega, beta, eta)
            
            #print("Differences between GF components for Klein terms t1 and t2 are", G_RR_t1 - G_RR_t2,  
            #      G_RL_t1 - G_RL_t2, G_LR_t1 - G_LR_t2, G_LL_t1 - G_LL_t2)

            toc =time.perf_counter()
            #print("time for summation: ", toc-tic)
            
            # print(summation_sector(f_B_R, f_B_R, energies_l1, energies_l2, omega, beta, eta))
            # write to GF array
            GF[l, :] += prefactor * np.array([G_RR_t1 + G_RR_t1, G_RL_t1 + G_RL_t2,
              G_LR_t1 + G_LR_t2, G_LL_t1 + G_LL_t2])

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
    
    omega= 0.
    beta = [1, 5, 10, 20]
    beta = [0.1]
    eta = [0.01]

    
    
    plot_spectrum(parameters, omega, beta, eta)


if __name__ == '__main__':
    main()




	