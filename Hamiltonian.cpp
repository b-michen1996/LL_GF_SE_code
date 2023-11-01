/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */

#include "Hamiltonian.h"
#include "multi_index_aux.h"

#include <cmath>
#include <array>
#include <iostream>
#include <omp.h>


M H_Luttinger_J(double v_F, double K, double U_m, double L, double a, double alpha, vector<vector<int>> HS_sector, int J){
    /* Calculates matrix representation of the total bosonic Luttinger Hamiltonian 
     * H0_B + J * H1_B on the subspace given by the Fock states / vectors of occupation 
     * numbers in HS_sector. */
    int dim_sector = HS_sector.size();
    
    // initialize matrix to store results
    M result = M::Zero(dim_sector,dim_sector);        
    
    complex<double> prefactor_H0 = v_F * 2 * M_PI / (L * K);    
    complex<double> prefactor_H1 = pow(2 * M_PI * alpha/ L, K - 1) * 1i * a * U_m/ L;    
    
    double norm_H0 = 0;
    double norm_H1 = 0;
    
    //omp_set_num_threads(14);
    #pragma omp parallel for
    for (int l_1 = 0; l_1 < dim_sector; l_1++){
        // set diagonal entry
        result(l_1, l_1) = prefactor_H0 * energy_H0(HS_sector[l_1]);
        norm_H0 += pow(abs(result(l_1, l_1)), 2); 
        // set off-diagonal entries
        for (int l_2 = l_1 + 1; l_2 < dim_sector; l_2++){
            vector<int> beta_1 = HS_sector[l_1];
            vector<int> beta_2 = HS_sector[l_2];
            result(l_1, l_2) = double(J) * prefactor_H1 * H1_B_matrix_element(K, beta_1, beta_2, L, a); 
            result(l_2, l_1) = -result(l_1, l_2);
            norm_H1 += 2 * pow(abs(result(l_1, l_2)), 2);
            }
        }
    
    cout << "norm H0 = " << sqrt(norm_H0) << ", norm H1 = " << sqrt(norm_H1) << "\n";
    
    return result;
};
        
        
double H1_B_matrix_element(double K, vector<int> beta_1, vector<int> beta_2, double L, 
        double a){
    /* Calculate matrix element of H1_B between two Fock states beta_1, beta_2*/
    
    // we can check immediately if the matrix element vanishes
    if ((abs_mi(beta_1) + abs_mi(beta_2)) % 2 == 0){
        return 0;
    }
    
    int dim = beta_1.size();
    double sqrt_K = sqrt(K);

    // We need an  array of lower and upper bounds for each component.
    vector<int> lower(dim);
    vector<int> upper(dim);
    
    for (int l= 0; l < dim; l++){  
        lower[l] = max(beta_1[l], beta_2[l]);
        upper[l] = beta_1[l] + beta_2[l];
    }
    
    // loop over all allowed values of alpha using auxilliary function
    vector<int> alpha = lower;
    double result = 0;
    double factor = 2 * M_PI * a / L;
    int p_c = int(dim / 2 + 0.1);
    
    while (alpha.size() > 0){
        // We need to perform sums of multi-indices element-wise
        vector<int> two_a_m_b1_m_b2(dim);
        vector<int> a_m_b1(dim);
        vector<int> a_m_b2(dim);
        vector<int> b1_p_b2_m_a(dim);
        for (int l= 0; l < dim; l++){ 
            a_m_b1[l] = alpha[l] - beta_1[l];
            a_m_b2[l] = alpha[l] - beta_2[l];
            two_a_m_b1_m_b2[l] = a_m_b1[l] + a_m_b2[l];
            b1_p_b2_m_a[l] = beta_1[l] + beta_2[l] - alpha[l];
        }

        // calculate term in sum
        result += power_sqrt_l_over_l_mi(two_a_m_b1_m_b2)  * (function_A(a_m_b1, factor, p_c) 
                - function_A(a_m_b2, factor, p_c)) * 2 * pow(sqrt_K, abs_mi(a_m_b1) 
                + abs_mi(a_m_b2)) / (factorial_mi(a_m_b1) 
                * factorial_mi(a_m_b2) * factorial_mi(b1_p_b2_m_a));   
        
        alpha = next_val(lower, upper, alpha);
    };    

    return sqrt(factorial_mi(beta_1)) * sqrt(factorial_mi(beta_2)) * result;
}


eigenstates_J_pm::eigenstates_J_pm(double in_v_F, double in_K, double in_U_m, double in_L, 
        double in_a, double in_alpha, int in_E_c, int in_p_c, int in_threads):v_F(in_v_F), 
        K(in_K), U_m(in_U_m), L(in_L), a(in_a), alpha(in_alpha), E_c(in_E_c), p_c(in_p_c), 
        threads(in_threads){
    /* Constructor for HS_Ec_pc struct containing all sectors of total Hilbert 
     * space that is truncated up to single-particle momentum p_c and energy E_c*/
    
    // Calculate all eigenstates and eigen-energies in each total momentum sector
    energies_ES_sector_wise_J_1.resize(2 * E_c + 1);
    energies_ES_sector_wise_J_m1.resize(2 * E_c + 1);
    
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (int l = 0; l < 2 * E_c + 1; l++){
        // momentum sector
        vector<vector<int>> HS_sector = HS_truncated.HS_tot[l];
        
        // Get matrix representation of bosonic Hamiltonian
        M H_Luttinger_J_1_block = H_Luttinger_J(v_F, K, U_m, L, a, alpha, HS_sector, 1);
        M H_Luttinger_J_m1_block = H_Luttinger_J(v_F, K, U_m, L, a, alpha, HS_sector, -1);
        
        cout << "total momentum sector " << l - E_c << ", norm of Hamiltonian "
                << H_Luttinger_J_1_block.norm() << "\n";
        
        // Calculate eigenstates and insert into global vector
        Eigen::SelfAdjointEigenSolver<M> temp_J1(H_Luttinger_J_1_block);    
        Eigen::SelfAdjointEigenSolver<M> temp_Jm1(H_Luttinger_J_m1_block);

        energies_ES_sector_wise_J_1[l] = temp_J1;
        energies_ES_sector_wise_J_m1[l]= temp_Jm1;
    }   
    
}