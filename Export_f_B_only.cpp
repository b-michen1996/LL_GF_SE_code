/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */
#include <iostream>
#include <omp.h>
#include <fstream>
#include <filesystem>

#include "GF_SE.h"
#include "multi_index_aux.h"

#include "Export_f_B_only.h"

void export_f_B(double v_F, double K, double U_m, double L, double a, double alpha, 
        int E_c, int p_c, int runNr, int threads){
    /* Calculate only the matrix elements of the bosonic factors f_B_r and 
     * export them together with the energies.*/
    
    // Create folder and save parameters
    string output_folder = "result_f_B_run_" + to_string(runNr);
    filesystem::create_directories(output_folder);
    ofstream parameters(output_folder + "/parameters.csv");    
    
    parameters << "v_F " << v_F << "\n" << "K " << K << "\n" << "alpha " << alpha << "\n" 
            << "U_m " << U_m << "\n" << "L " << L << "\n" 
            << "a " << a << "\n" << "E_c " << E_c << "\n" << "p_c " << p_c << "\n"  
            << "runNr " << runNr << "\n" << "threads " << threads;
    parameters.close();
    
    // Calculate eigenstates
    eigenstates_J_pm ES(v_F, K, U_m, L, a, alpha, E_c, p_c, threads);
    
    // export energies
    string energy_folder = output_folder + "/energies";
    filesystem::create_directories(energy_folder);

    for (int l = 0; l < 2 * E_c + 1; l++){        
        // get dimension of current momentum sector
        int dim_sector = ES.HS_truncated.HS_tot[l].size();
        // get current sector
        Eigen::SelfAdjointEigenSolver<M> solver_curr_J_1 = ES.energies_ES_sector_wise_J_1[l]; 
        Eigen::SelfAdjointEigenSolver<M> solver_curr_J_m1 = ES.energies_ES_sector_wise_J_m1[l]; 
        
        // create files 
        ofstream energies_J_1(energy_folder + "/energies_J_1_" + to_string(l) + ".csv"); 
        ofstream energies_J_m1(energy_folder + "/energies_J_m1_" + to_string(l) + ".csv"); 
        
        for (int j = 0; j < dim_sector; j++){
            double E_J_1_j = solver_curr_J_1.eigenvalues()(j);
            double E_J_m1_j = solver_curr_J_m1.eigenvalues()(j);           
            energies_J_1 << E_J_1_j << "\n";
            energies_J_m1 << E_J_m1_j << "\n";
        } 
        energies_J_1.close();
        energies_J_m1.close();
    }
    
    
    
    // export f_B  
    vector<vector<int>> HS_sector_l1 = ES.HS_truncated.HS_tot[4];
    vector<vector<int>> HS_sector_l2 = ES.HS_truncated.HS_tot[3];  
        
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (int l = 0; l < 2 * p_c + 1; l++ ){
        int k = l - p_c;
        
        // create subfolder to store f_B(k)
        string output_folder_f_B_R_curr = output_folder + "/f_B_R_" + to_string(l);
        string output_folder_f_B_L_curr = output_folder + "/f_B_L_" + to_string(l);
        filesystem::create_directories(output_folder_f_B_R_curr);
        filesystem::create_directories(output_folder_f_B_L_curr);
        
        
        for (int l_2 = max(0, k); l_2 < 2 * E_c + 1; l_2 ++){
            int l_1 = l_2 - k;
            if (l_1 > 2 * E_c){
                break;
            }
            // Get eigenstates and energies of both momentum and J sectors 
            Eigen::SelfAdjointEigenSolver<M> Eigen_J_1_l1 = ES.energies_ES_sector_wise_J_1[l_1];
            Eigen::SelfAdjointEigenSolver<M> Eigen_J_1_l2 = ES.energies_ES_sector_wise_J_1[l_2];
            Eigen::SelfAdjointEigenSolver<M> Eigen_J_m1_l1 = ES.energies_ES_sector_wise_J_m1[l_1];
            Eigen::SelfAdjointEigenSolver<M> Eigen_J_m1_l2 = ES.energies_ES_sector_wise_J_m1[l_2];

            // get matrix rep of f_B_r(k) in fock basis
            vector<vector<int>> HS_sector_l1 = ES.HS_truncated.HS_tot[l_1];
            vector<vector<int>> HS_sector_l2 = ES.HS_truncated.HS_tot[l_2];        
                
            M f_B_R_curr_Fock_base = f_B_r(1, K, HS_sector_l1, HS_sector_l2);
            M f_B_L_curr_Fock_base = f_B_r(-1, K, HS_sector_l1, HS_sector_l2);  
            
            // do matrix multiplication to obtain all scalar products
            M f_B_R_curr_ES_base = (Eigen_J_1_l1.eigenvectors().adjoint()) 
                    * (f_B_R_curr_Fock_base * Eigen_J_1_l2.eigenvectors());
            M f_B_L_curr_ES_base = (Eigen_J_1_l1.eigenvectors().adjoint()) 
                    * (f_B_L_curr_Fock_base * Eigen_J_1_l2.eigenvectors());
            
            // write to file
            string filename_f_B_R= output_folder_f_B_R_curr + "/f_B_R_" + to_string(l_1) 
            + "_" + to_string(l_2) + ".csv"; 
            string filename_f_B_L = output_folder_f_B_L_curr + "/f_B_L_" + to_string(l_1) 
            + "_" + to_string(l_2) + ".csv";
            
            cout << "l_1 = " << l_1 << ", l_2 = " << l_2 << ", norm f_B_R_curr_Fock_base = " 
                    << f_B_R_curr_Fock_base.norm() << ", norm f_B_L_curr_Fock_base = " 
                    << f_B_L_curr_Fock_base.norm() << "\n";
            
            write_f_B_to_file(f_B_R_curr_ES_base, filename_f_B_R);
            write_f_B_to_file(f_B_L_curr_ES_base, filename_f_B_L);
        }
    }
}


void write_f_B_to_file(M f, string filename){
    /* Write matrix to file in a format that numpy can read.*/
    
    // create file to store result
    ofstream file(filename); 
  
    int rows = f.rows();
    int cols = f.cols();
 
    for (int j_1 = 0; j_1 < rows;  j_1++){  
        for (int j_2 = 0; j_2 < cols;  j_2++){  
            double M_E_real = f(j_1, j_2).real();
            double M_E_imag = f(j_1, j_2).imag();
            
            if (M_E_imag > 0){
                file << M_E_real << "+" << abs(M_E_imag) << "j ";
            }
            else{
                file << M_E_real << "-" << abs(M_E_imag) << "j ";
            }             
        }
        file << "\n";
    }
    file.close();
}

