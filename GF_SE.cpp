/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */


#include <iostream>
#include <omp.h>
#include <fstream>
#include <filesystem>

#include "GF_SE.h"

void GF_SE_J_explicit(double v_F, double K, double U_m, double L, double a, double alpha, int E_c, int p_c, double omega,
        int runNr){
    /* Calculate retarded GF and associated self-energy at omega and save to file*/
    
    // Create folder and save parameters
    string output_folder = "result_run_" + to_string(runNr);
    filesystem::create_directories(output_folder);
    ofstream result(output_folder + "/result.csv");
    ofstream parameters(output_folder + "/parameters.csv");    
    
    parameters << "v_F " << v_F << "\n" << "K " << K << "\n" << "alpha " << alpha << "\n" 
            << "U_m " << U_m << "\n" << "L " << L << "\n" 
            << "a " << a << "\n" << "E_c " << E_c << "\n" << "p_c " << p_c << "\n"  
            << "runNr " << runNr << "\n" ;
    parameters.close();
    
    return;
}


eigenstates_J_pm::eigenstates_J_pm(double in_v_F, double in_K, double in_U_m, double in_L, 
        double in_a, double in_alpha, int in_E_c, int in_p_c):v_F(in_v_F), K(in_K), U_m(in_U_m), L(in_L), 
        a(in_a), alpha(in_alpha), E_c(in_E_c), p_c(in_p_c){
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
        M H_Luttinger_J_1_block = H_Luttinger_J(v_F, K, U_m, L, a, alpha, HS_sector, 1);
        M H_Luttinger_J_m1_block = H_Luttinger_J(v_F, K, U_m, L, a, alpha, HS_sector, -1);
        
        // Calculate eigenstates and insert into global vector
        Eigen::SelfAdjointEigenSolver<M> temp_J1(H_Luttinger_J_1_block);    
        Eigen::SelfAdjointEigenSolver<M> temp_Jm1(H_Luttinger_J_m1_block);
        cout << "Set of ES and energies calculated \n";
        energies_ES_sector_wise_J_1.insert(energies_ES_sector_wise_J_1.end(),temp_J1);
        energies_ES_sector_wise_J_m1.insert(energies_ES_sector_wise_J_m1.end(),temp_Jm1);        
    }
}


double f_B_r_matrix_element(int r, vector<int> beta_1, vector<int> beta_2){
    /* Calculates the matrix element of the bosonic factor f_B_r(k) (the dependence 
     * on the external momentum k is only in the total momentum sectors between 
     * which we evaluate f_B_r(k)).*/
    double result = 0;
    
    return result;
}