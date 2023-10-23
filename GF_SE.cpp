/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */


#include <iostream>
#include <omp.h>

#include "GF_SE.h"


eigenstates_J_pm::eigenstates_J_pm(double in_u, double in_U_m, double in_L, 
        double in_a, int in_E_c, int in_p_c):u(in_u), U_m(in_U_m), L(in_L), 
        a(in_a), E_c(in_E_c), p_c(in_p_c){
    /* Constructor for HS_Ec_pc struct containing all sectors of total Hilbert 
     * space that is truncated up to single-particle momentum p_c and energy E_c*/
    
    // generate truncated Hilbert space
    //HS_Ec_pc HS_truncated(E_c, p_c);
    
    // Calculate all eigenstates and eigen-energies in each total momentum sector
    
    omp_set_num_threads(14);
    #pragma omp parallel for
    for (int l = 0; l < 2 * E_c + 1; l++){
        // momentum sector
        vector<vector<int>> HS_sector = HS_truncated.HS_tot[l];
        
        // Get matrix representation of bosonic Hamiltonian
        M H_Luttinger_J_1_block = H_Luttinger_J(u, U_m, L, a, HS_sector, 1);
        M H_Luttinger_J_m1_block = H_Luttinger_J(u, U_m, L, a, HS_sector, -1);
        
        // Calculate eigenstates and insert into global vector
        Eigen::SelfAdjointEigenSolver<M> temp_J1(H_Luttinger_J_1_block);    
        Eigen::SelfAdjointEigenSolver<M> temp_Jm1(H_Luttinger_J_m1_block);
        cout << "Set of ES and energies calculated \n";
        energies_ES_sector_wise_J_1.insert(energies_ES_sector_wise_J_1.end(),temp_J1);
        energies_ES_sector_wise_J_m1.insert(energies_ES_sector_wise_J_m1.end(),temp_Jm1);        
    }
}