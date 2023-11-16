/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/main.cc to edit this template
 */

/* 
 * File:   main.cpp
 * Author: bm
 *
 * Created on October 18, 2023, 3:59 PM
 */

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <chrono>

#include "hilbert_space.h"
#include "Hamiltonian.h"
#include "multi_index_aux.h"
#include "GF_SE.h"
#include "Export_f_B_only.h"
#include "export_b_q_only.h"

using namespace std;
using namespace std::chrono;

typedef Eigen::MatrixXcd M;
/*
 * 
 */
int main(int argc, char** argv){   
    // some parameters
    double v_F = 0.5;
    double K = 1;
    double L = 100;
    double a = 1;        
    double U_m = 0;
    int E_c = 12;
    
    double gamma = 0.001;
    double eta = 0.01;
    
    double omega = 0.00;
    double beta = 10;
    
    int runNr = 2;
    int threads = 15;
    
    int p_c = 12;
    
    /*
    GF_SE_J_explicit(v_F, K, U_m, L, a, gamma, eta, E_c, p_c, omega, 
            beta, runNr, threads);
    */
    
    /*
    export_f_B(v_F, K, U_m, L, a, gamma, E_c, p_c,
            runNr, threads);
    */
    /*
    export_f_B_Klein_explicit(v_F, K, U_m, L, a, gamma, E_c, p_c,
            runNr, threads);
    */
  
    
    export_b_q_only(v_F, K, U_m, L, a, gamma, E_c, p_c,
            runNr, threads);
    
    /*
    HS_Ec_pc HS_test{E_c, p_c};
    
    for (int j = 0; j < 2 * E_c + 1; j++) {
        vector<vector<int>> HS_sector = HS_test.HS_tot[j]; 
        cout << "j = " << j << ", momentum = " <<  j - E_c << ", states :\n";
        
        for (int l = 0; l < HS_sector.size(); l++){
            vector<int> alpha_curr = HS_sector[l];
            for (int x : alpha_curr){
                cout << " " << x;
            }
            cout << "\n";
        }
        
                
    }
    */
    return 0;
}


/*
vector<int> alpha_test = {0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0};
    
    cout << "P(alpha) = " << momentum_mi(alpha_test) << "\n"; 
    cout << "energy(alpha) = " << energy_H0_reg(10 * L, alpha_test, gamma) << "\n"; 
    
    
    // Calculate eigenstates
    eigenstates_J_pm ES(v_F, K, U_m, L, a, gamma, E_c, p_c, threads);
    
    int l1 = 4;
    int l2 = 8;
    
    // Get eigensolver for l
    Eigen::SelfAdjointEigenSolver<M> solver_curr_J_1_l1 = ES.energies_ES_sector_wise_J_1[l1]; 
    Eigen::SelfAdjointEigenSolver<M> solver_curr_J_1_l2 = ES.energies_ES_sector_wise_J_1[l2]; 
    
    // Hamiltonian matrix
    vector<vector<int>> HS_sector = ES.HS_truncated.HS_tot[l1];    
    M H_Luttinger_J_1_block = H_Luttinger_J(v_F, K, U_m, L, a, gamma, HS_sector, 1);
    
    // extract some eigenvector and eigenvalue for testing
    int k = 3;
    
    Eigen::VectorXcd eigenvector_test = solver_curr_J_1_l1.eigenvectors().col(k);    
    double energy_test = solver_curr_J_1_l1.eigenvalues()(k);
    
    //cout << "energy: " << energy_test << "\n";
    //cout << "eigenvector " <<  "\n" << eigenvector_test <<  "\n";
    //cout << "test eigenvalue: " << eigenvector_test.dot(H_Luttinger_J_1_block * eigenvector_test) << "\n";
    //cout << "test overlap" << pow((solver_curr_J_1_l1.eigenvectors().adjoint() 
            //* solver_curr_J_1_l2.eigenvectors()).norm(), 2) << "\n";
    cout << "energies l1 = " << l1 << "\n" << solver_curr_J_1_l1.eigenvalues() <<  "\n";
    cout << "energies l2 = " << l2 << "\n" << solver_curr_J_1_l2.eigenvalues() <<  "\n";
    
    double prefactor = gamma * 2 * M_PI / L;
    
    for (int j = 0; j < HS_sector.size(); j++){
        vector<int> alpha_curr = HS_sector[j];
        cout << "-------------- \n";
        cout << "j = " << j << ", energy = " << energy_H0(alpha_curr) 
                << ", energy_reg = " << energy_H0_reg(L, alpha_curr, gamma) << "\n";
        cout << "alpha = ";
        
        for (int x : alpha_curr){
            cout << " " << x;
        }
        cout << "\n P(alpha) = " << momentum_mi(alpha_curr);
        
        Eigen::VectorXcd eigenvector_curr = solver_curr_J_1_l2.eigenvectors().col(j);
                
        cout << "\n Eigenstate norm " << eigenvector_curr.norm() <<"\n";
        
        for (int m = 0; m < eigenvector_curr.size() ; m++){
            cout << eigenvector_curr[m] << " ";
        }
        
        cout << "\n";
    }
    */

    
    



