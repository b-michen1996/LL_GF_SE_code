/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */
#include <iostream>
#include <omp.h>
#include <fstream>
#include <filesystem>

#include "export_b_q_only.h"
#include "Export_f_B_only.h"
#include "multi_index_aux.h"

void export_b_q_only(double v_F, double K, double U_m, double L, double a, double gamma, 
        int E_c, int p_c, int runNr, int threads){
    /* Calculate only the matrix elements of the creators and annihilators b_q, b_q^\dagger 
     * and export them together with the energies.*/
    
    // Create folder and save parameters
    string output_folder = "result_b_q_run_" + to_string(runNr);
    filesystem::create_directories(output_folder);
    ofstream parameters(output_folder + "/parameters.csv");    
    
    parameters << "v_F " << v_F << "\n" << "K " << K << "\n" << "gamma " << gamma << "\n" 
            << "U_m " << U_m << "\n" << "L " << L << "\n" 
            << "a " << a << "\n" << "E_c " << E_c << "\n" << "p_c " << p_c << "\n"  
            << "runNr " << runNr << "\n" << "threads " << threads;
    parameters.close();
    
    // Calculate eigenstates
    eigenstates_J_pm ES(v_F, K, U_m, L, a, gamma, E_c, p_c, threads);
    
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
            double E_J_1_j = solver_curr_J_1.eigenvalues()[j];
            double E_J_m1_j = solver_curr_J_m1.eigenvalues()[j];           
            energies_J_1 << E_J_1_j << "\n";
            energies_J_m1 << E_J_m1_j << "\n";
        } 
        energies_J_1.close();
        energies_J_m1.close();
    }
      
    // export b_q          
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (int l = 0; l < 2 * p_c; l++ ){
        int k = momentum(l, p_c);
        
        // create subfolder to store b_q, b_q^dagger
        string output_folder_b_q_curr = output_folder + "/b_" + to_string(l);

        filesystem::create_directories(output_folder_b_q_curr);        
        
        for (int l_2 = max(0, k); l_2 < 2 * E_c + 1; l_2 ++){
            int l_1 = l_2 - k;
            if (l_1 > 2 * E_c){
                break;
            }
            // Get eigenstates and energies of both momentum sectors 
            Eigen::SelfAdjointEigenSolver<M> Eigen_l1 = ES.energies_ES_sector_wise_J_1[l_1];
            Eigen::SelfAdjointEigenSolver<M> Eigen_l2 = ES.energies_ES_sector_wise_J_1[l_2];

            // get matrix rep of b_q in Fock basis
            vector<vector<int>> HS_sector_l1 = ES.HS_truncated.HS_tot[l_1];
            vector<vector<int>> HS_sector_l2 = ES.HS_truncated.HS_tot[l_2];        
            
            M b_q_curr_Fock_base = b_q_Fock(l, HS_sector_l1, HS_sector_l2);
            
            // do matrix multiplication to obtain all scalar products
            M b_q_curr_ES_base = ((Eigen_l1.eigenvectors()).adjoint()) 
                    * (b_q_curr_Fock_base * Eigen_l2.eigenvectors());
                        
            //cout << Eigen_l1.eigenvectors();
            // write to file
            string filename_b_q = output_folder_b_q_curr + "/b_" + to_string(l_1) 
            + "_" + to_string(l_2) + ".csv"; 
            
            write_f_B_to_file(b_q_curr_ES_base, filename_b_q);
        }
    }
}


M b_q_Fock(int l, vector<vector<int>> HS_sector_1, vector<vector<int>> HS_sector_2){
    /* Calculates Fock basis matrix representation of b_q between two total momentum sectors.*/
    int dim_sector_1 = HS_sector_1.size();
    int dim_sector_2 = HS_sector_2.size();
    int dim = HS_sector_1[0].size();
    
    // initialize matrix to store results
    M result = M::Zero(dim_sector_1, dim_sector_2);        

    for (int l_1 = 0; l_1 < dim_sector_1; l_1++){        
        for (int l_2 = 0; l_2 < dim_sector_2; l_2++){
            vector<int> beta_1 = HS_sector_1[l_1];
            vector<int> beta_2 = HS_sector_2[l_2];
            
            double temp = 1;
            // loop over all entries and check if they are the same apart 
            // from the one where b_q acts
            for (int j = 0; j < dim; j++){
                // skip if j = l
                if (j == l){
                    continue;
                }
                
                // if the other entries do not agree, set result to zer
                
                if (!(beta_2[j] == beta_1[j])){
                    temp *= 0.; 
                }
            }
            
            if (beta_2[l] - 1 == beta_1[l]){
                temp *=  sqrt(beta_2[l]);
            }else {
                temp *= 0;
            }   
            
            // set matrix element
            result(l_1, l_2) = temp;
            }
        }
    return result;
};

