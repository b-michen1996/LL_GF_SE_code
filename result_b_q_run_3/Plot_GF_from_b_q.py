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


def plot_GF_phi(parameters, omega_1, omega_2, beta, eta):
    """Plot syndromes data"""
    # get some parameters
    v_F = parameters[0,1]
    K = parameters[1,1]
    alpha = parameters[2,1]
    U_m = parameters[3,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
    
    title = (r"$v_F$ = " + str(v_F) +  r", K = " + str(K) +  r", $U_m$ = " 
             + str(U_m) + r", $\gamma$ = " + str(alpha) + r", L = " + str(L) 
    + r", $E_c$ = " + str(E_c) + r", $p_c$ =" + str(p_c)) 
    
    # plot
    fig = plt.figure()
    a1 =  fig.add_subplot(1,1,1)
    a1.set_title(title , fontsize=30)
    a1.set_xlabel("k", fontsize=30)
    a1.set_ylabel("E", fontsize=30)

    # calculate GF
    GF = calc_GF_phi(parameters, omega_1, omega_2, beta)
    
    plot_data = 1j * np.zeros(2* p_c + 1)
    plot_data[: p_c] = GF[0: p_c, 0]
    plot_data[p_c] = 1j * 0
    plot_data[p_c + 1:] = GF[p_c :, 0]
    
    plot_data_R = 1j * np.zeros(2* p_c + 1)
    plot_data_R[: p_c] = GF[0: p_c, 1]
    plot_data_R[p_c] = 1j * 0
    plot_data_R[p_c + 1:] = GF[p_c :, 1]
    
    plot_data_analytical = GF_phi_analytically(parameters, omega_1, beta)
    print(plot_data_analytical )

    # calculate momenta
    k_vals = np.linspace(-p_c, p_c, 2 * p_c + 1) * 2 * np.pi / L
    
    a1.plot(k_vals, plot_data.real, color = "r", linestyle = "-", label = "Numerics Matsubara") 
    a1.plot(k_vals, plot_data.imag, color = "r",  linestyle = "--")
    
    a1.plot(k_vals, plot_data_R.real, color = "g", linestyle = "-", marker = "x",
            label = "Numerics retarded") 
    a1.plot(k_vals, plot_data_R.imag, color = "g",  linestyle = "--")
    
    a1.plot(k_vals, plot_data_analytical.real, color = "b", linestyle = "-", label = "Analytical") 
    a1.plot(k_vals, plot_data_analytical.imag, color = "b", linestyle = "--") 

    #a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    a1.legend(loc='upper center', prop={'size': 12})
    # change color of legend to black
    leg = a1.get_legend()
    
    plt.show()
    
    return
  
    
def plot_spectrum_GF_bq(parameters, omega, beta, eta):
    """Plot spectrum of self-energy obtained from b_q correlator"""
    # get some parameters
    v_F = parameters[0,1]
    K = parameters[1,1]
    alpha = parameters[2,1]
    U_m = parameters[3,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
    

    title = (r"$v_F$ = " + str(v_F) +  r", K = " + str(K) +  r", $U_m$ = " + str(U_m) + r", $\gamma$ = " + str(alpha)
    + r", L = " + str(L)) + r", $E_c$ = " + str(E_c) + r", $p_c$ =" + str(p_c) 
    # plot
    fig = plt.figure()
    a1 =  fig.add_subplot(1,1,1)
    a1.set_title(title , fontsize=30)
    a1.set_xlabel("k", fontsize=30)
    a1.set_ylabel("E", fontsize=30)
    
    # calculate GF
    GF = calc_GF_from_bq(parameters, omega, beta, eta)

    # calculate spectrum    
    spectrum = calc_spectrum_from_GF(GF, omega, p_c)

    # calculate momenta
    k_vals = np.linspace(-p_c, p_c, 2 * p_c + 1) * 2 * np.pi / L
    
    # calculate expectation
    spectrum_expected = np.abs(k_vals) * (v_F / (2 * K))
    
    
    a1.plot(k_vals, spectrum.real, color = "r", linestyle = "-",
            label = r"$\omega$ = " + str(omega) + r", $\beta$ = " + str(beta) + r", $\eta$ = " + str(eta) ) 
    a1.plot(k_vals, spectrum.imag, color = "r",  linestyle = "--")
    
    a1.plot(k_vals, spectrum_expected, color = "b", linestyle = "-",
            label = "Expected spectrum" ) 

    #a1.set_xlim([-np.pi, np.pi])
    #a1.set_ylim([-2.5, 2.5])
    a1.tick_params(direction='out', length=6, width=2,  labelsize = 'large')
    a1.legend(loc='upper center', prop={'size': 12})
    # change color of legend to black
    leg = a1.get_legend()
    
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


            
def calc_GF_phi(parameters, omega_1, omega_2, beta):
    """ Read in matrix-representation of bosonic factors form file and calculate
        GF of real bosonic field phi for given omega, beta, eta."""
    # get some parameters
    K = parameters[1,1]
    alpha = parameters[2,1]
    L = parameters[4,1]
    E_c = int(parameters[6,1])
    p_c = int(parameters[7,1]) 
        
    GF = 1j * np.zeros((int(2 * p_c), 2))
    
    # calculate partition sum and global prefactor
    Z = 0
    for l in range(0, 2 * E_c + 1): 
        # get energies of sector
        filename_energies_l = "energies/energies_J_1_" + str(l) + ".csv" 
        energies_l = np.genfromtxt(filename_energies_l, dtype=float, delimiter=' ')
                
        # add to sum 
        Z += np.sum(np.exp(- beta * energies_l))
        
    print("partition sum Z = ", Z)
    prefactor = K * np.pi / (2 * L * Z)


    # loop over all momenta from -p_c to p_c
    for l in range(0, 2 * p_c):
        q = l - p_c + 1
        if (l < p_c):
            q = l - p_c
        
        # loop over all momentum sectors with non-vanishing matrix elements of b_q        
        for l2 in range(max(0, -q) , 2 * E_c + 1):
            l1 = l2 + q              
            if l1 >  2 * E_c:
                break
            tic =time.perf_counter()
            
            # get index of negative momentum
            ml = 2 * p_c - l - 1
            
            # get matrix representation of b_q in eigenbasis
            filename_b_ql_l2_l1 = "b_" + str(l) + "/b_" + str(l2) + "_" + str(l1) + ".csv"
            filename_b_mql_l1_l2 = "b_" + str(ml) + "/b_" + str(l1) + "_" + str(l2) + ".csv"
            
            b_ql_l2_l1 = np.genfromtxt(filename_b_ql_l2_l1, dtype=complex, delimiter=' ')
            b_mql_l1_l2 = np.genfromtxt(filename_b_mql_l1_l2, dtype=complex, delimiter=' ')
            
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
            G_sector = summation_sector_phi(b_ql_l2_l1, b_mql_l1_l2, energies_l1, energies_l2, 
                                       omega_1, omega_2, beta)

            G_R_sector = summation_R_sector_phi(b_ql_l2_l1, b_mql_l1_l2, energies_l1, energies_l2, 
                                       omega_1, beta)


            toc =time.perf_counter()
            #print("time for summation: ", toc-tic)
            
            # print(summation_sector(f_B_R, f_B_R, energies_l1, energies_l2, omega, beta, eta))
            # write to GF array
            GF[l, 0] += (prefactor / abs(q * 2 * np.pi / L)) * G_sector
            GF[l, 1] += (prefactor / abs(q * 2 * np.pi / L)) * G_R_sector

    return GF


def GF_phi_analytically(parameters, omega, beta):
    """Analytical expression for GF"""
    v_F = parameters[0,1]
    K = parameters[1,1]
    L = parameters[4,1]
    p_c = int(parameters[7,1])
    
    u = v_F / K
    
    result = 1j * np.zeros(2 * p_c + 1)
    
    for j in range(0, p_c):
        q_j = (j - p_c) * 2 * np.pi / L
        result[j] = (np.pi * K * L * beta) / (omega**2 / u + u * q_j ** 2)
    
    for j in range(p_c + 1, 2 * p_c + 1):
        q_j = (j - p_c) * 2 * np.pi / L
        result[j] = (np.pi * K * L * beta) / (omega**2 / u + u * q_j**2)
    
    return result


def summation_sector_phi(b_ql_l2_l1, b_mql_l1_l2, energies_l1, energies_l2, 
                        omega_1, omega_2, beta):
    """"Carry out the summation from the Kaellen-Lehmann representation for 
    the given parameters"""
    numerator = (np.exp(-beta * energies_l1) * (np.ones(energies_l1.shape) 
    - np.exp(beta * (-1j * omega_1 * np.ones(energies_l1.shape) + energies_l1 
    - energies_l2))) * (np.ones(energies_l1.shape) 
    - np.exp(beta * (1j * omega_2 * np.ones(energies_l1.shape) + energies_l2 
    - energies_l1))))
    
    denominator = (-1j *omega_1 * np.ones(energies_l1.shape) + energies_l1 
    - energies_l2) * (1j *omega_2 * np.ones(energies_l1.shape) + energies_l2 
    - energies_l1)
    
    factor_matrix = np.divide(numerator, denominator)
        
    temp = (b_mql_l1_l2 * np.conj(b_mql_l1_l2) + b_mql_l1_l2 * np.transpose(b_ql_l2_l1)
    + np.conj(b_ql_l2_l1.T)* np.conj(b_mql_l1_l2) 
    + np.conj(b_ql_l2_l1.T) * np.transpose(b_ql_l2_l1))
    
    result = np.sum(temp * factor_matrix)
    
    return result


def summation_R_sector_phi(b_ql_l2_l1, b_mql_l1_l2, energies_l1, energies_l2, 
                        omega, beta):
    """"Carry out the summation from the Kaellen-Lehmann representation for 
    the given parameters"""
    numerator = np.exp(-beta * energies_l1)
    
    denominator = (1j * omega * np.ones(energies_l1.shape) + energies_l1 
    - energies_l2)
    
    factor_matrix = np.divide(numerator, denominator)
        
    temp = (b_mql_l1_l2 * np.conj(b_mql_l1_l2) + b_mql_l1_l2 * np.transpose(b_ql_l2_l1)
    + np.conj(b_ql_l2_l1.T)* np.conj(b_mql_l1_l2) 
    + np.conj(b_ql_l2_l1.T) * np.transpose(b_ql_l2_l1))
    
    result = np.sum(temp * factor_matrix)
    
    return result


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
    
    omega = 0.
    beta = [1, 5, 10, 20]
    beta = 0.5
    eta = 0.01

    p_c = 8
    
    n1 = 0
    n2 = 0
    
    omega_1 = -1j * (n1 * 2 * np.pi / beta + 1j * eta)
    omega_2 = -1j * (n1 * 2 * np.pi / beta + 1j * eta)
    
    #omega_1 = (n1 * 2 * np.pi / beta + 1j * eta)
    #omega_2 = (n1 * 2 * np.pi / beta + 1j * eta)
    
    #GF_bosonic = calc_GF_from_b_q(parameters, omega_1, omega_2, beta)
    
    #plot_GF_phi(parameters, omega_1, omega_2, beta, eta)
    
    #plot_spectrum(parameters, omega, beta, eta)
    plot_spectrum_GF_bq(parameters, omega, beta, eta)


if __name__ == '__main__':
    main()




	