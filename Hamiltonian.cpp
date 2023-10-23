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


M H1_B(double L, double a, vector<vector<int>> HS_sector){
    /* Calculates matrix representation of bosonic Hamiltonian H1_B on the
     * subspace given by the Fock states / vectors of occupation numbers in 
     * HS_sector. */
    int dim_sector = HS_sector.size();
    
    // initialize matrix to store results
    M result = M::Zero(dim_sector,dim_sector);        
    
    complex<double> prefactor = 1i * a / (2 * M_PI * L);    
    
    omp_set_num_threads(14);
    #pragma omp parallel for
    for (int l_1 = 0; l_1 < dim_sector; l_1++){
        for (int l_2 = l_1 + 1; l_2 < dim_sector; l_2++){
            vector<int> beta_1 = HS_sector[l_1];
            vector<int> beta_2 = HS_sector[l_2];
            result(l_1, l_2) = prefactor * H1_B_matrix_element(beta_1, beta_2, L, a); 
            result(l_2, l_1) = -prefactor * result(l_1, l_2);
            }
        }
    
    
    return result;
};


M H_Luttinger_J(double u, double U_m, double L, double a, vector<vector<int>> HS_sector, int J){
    /* Calculates matrix representation of the total bosonic Luttinger Hamiltonian 
     * H0_B + J * H1_B on the subspace given by the Fock states / vectors of occupation 
     * numbers in HS_sector. */
    int dim_sector = HS_sector.size();
    
    // initialize matrix to store results
    M result = M::Zero(dim_sector,dim_sector);        
    
    complex<double> prefactor_H0 = u * 2 * M_PI / L;    
    complex<double> prefactor_H1 = 1i * a * U_m/ L;    
    
    //omp_set_num_threads(14);
    #pragma omp parallel for
    for (int l_1 = 0; l_1 < dim_sector; l_1++){
        // set diagonal entry
        result(l_1, l_1) = u * energy_H0(HS_sector[l_1]);
        
        // set off-diagonal entries
        for (int l_2 = l_1 + 1; l_2 < dim_sector; l_2++){
            vector<int> beta_1 = HS_sector[l_1];
            vector<int> beta_2 = HS_sector[l_2];
            result(l_1, l_2) = double(J) * prefactor_H1 * H1_B_matrix_element(beta_1, beta_2, L, a); 
            result(l_2, l_1) = -double(J) * prefactor_H1 * result(l_1, l_2);
            }
        }
    return result;
};
        
        
double H1_B_matrix_element(vector<int> beta_1, vector<int> beta_2, double L, 
        double a){
    /* Calculate matrix element of H1_B between two Fock states beta_1, beta_2*/
    
    // we can check immediately if the matrix element vanishes
    if ((abs_mi(beta_1) + abs_mi(beta_2)) / 2 == 0){
        return 0;
    }
    
    int dim = beta_1.size();

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
    int p_c = int(dim / 2);
    
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
        result += power_sqrt_l_over_l_mi(two_a_m_b1_m_b2) * (function_A(a_m_b1, factor, p_c) 
                - function_A(a_m_b2, factor, p_c)) * (1 - pow(-1, abs_mi(a_m_b1) 
                + abs_mi(a_m_b2))) / (factorial_mi(a_m_b1) 
                * factorial_mi(a_m_b2) * factorial_mi(b1_p_b2_m_a));   
        
        /*
        cout << "\n alpha = ";
            for (int x : alpha){
            cout << x << " ";
        }
        
        cout << "\n result = " << result << "\n";
        cout << " power_sqrt_l_over_l_mi(two_a_m_b1_m_b2) = " << power_sqrt_l_over_l_mi(two_a_m_b1_m_b2) << "\n";
        cout << " test = " << abs_mi(a_m_b1) + abs_mi(a_m_b2) << "\n";
        */
        alpha = next_val(lower, upper, alpha);
    };    
    
    /*
    cout << "\n factorial b1 = " << (factorial_mi(beta_1)) << "\n";
    cout << "\n factorial b2 = " << (factorial_mi(beta_2)) << "\n";
    cout << "\n sqrt (b1! b2!) = " << sqrt(factorial_mi(beta_1) * factorial_mi(beta_2)) << "\n";
    */
    
    return sqrt(factorial_mi(beta_1)) * sqrt(factorial_mi(beta_2)) * result;
}
