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

void GF_SE_J_explicit(double v_F, double K, double U_m, double L, double a, double alpha, 
        double eta, int E_c, int p_c, double omega, double beta, int runNr, int threads){
    /* Calculate retarded GF and associated self-energy at omega and save to file*/
    
    // Create folder and save parameters
    string output_folder = "result_run_" + to_string(runNr);
    filesystem::create_directories(output_folder);
    ofstream result(output_folder + "/result.csv");
    ofstream parameters(output_folder + "/parameters.csv");    
    
    parameters << "v_F " << v_F << "\n" << "K " << K << "\n" << "alpha " << alpha << "\n" 
            << "eta " << eta << "\n" << "U_m " << U_m << "\n" << "L " << L << "\n" 
            << "a " << a << "\n" << "E_c " << E_c << "\n" << "p_c " << p_c << "\n"  
            << "beta " << beta << "\n" << "runNr " << runNr << "\n" << "threads " << threads;
    parameters.close();
    
    // Calculate eigenstates
    eigenstates_J_pm ES(v_F, K, U_m, L, a, alpha, E_c, p_c, threads);
        
    // Calculate partition function
    double Z = 0;
    
    for (int l = 0; l < 2 * p_c + 1; l++){
        double Z_sector = 0;
        // get dimension of current momentum sector
        int dim_sector = ES.HS_truncated.HS_tot[l].size();
        // get current sector
        Eigen::SelfAdjointEigenSolver<M> solver_curr = ES.energies_ES_sector_wise_J_1[l];        
        for (int j = 0; j < dim_sector; j++){
            double E_j = solver_curr.eigenvalues().coeffRef(j);
            Z_sector += exp(- beta * E_j);
        }
        Z += Z_sector;        
    }
    
    // calculate elements of retarded GF for each k from -p_c to p_c 
    vector<complex<double>> G_RR;
    vector<complex<double>> G_RL;
    vector<complex<double>> G_LR;
    vector<complex<double>> G_LL;
    
    for (int l = 0; l < 2 * p_c + 1; l++ ){
        //int k = l - p_c;
    } 
    
    
    return;
}


eigenstates_J_pm::eigenstates_J_pm(double in_v_F, double in_K, double in_U_m, double in_L, 
        double in_a, double in_alpha, int in_E_c, int in_p_c, int in_threads):v_F(in_v_F), K(in_K), U_m(in_U_m), L(in_L), 
        a(in_a), alpha(in_alpha), E_c(in_E_c), p_c(in_p_c), threads(in_threads){
    /* Constructor for HS_Ec_pc struct containing all sectors of total Hilbert 
     * space that is truncated up to single-particle momentum p_c and energy E_c*/
    
    // generate truncated Hilbert space
    //HS_Ec_pc HS_truncated(E_c, p_c);
    
    // Calculate all eigenstates and eigen-energies in each total momentum sector
    
    omp_set_num_threads(threads);
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
        energies_ES_sector_wise_J_1.insert(energies_ES_sector_wise_J_1.end(),temp_J1);
        energies_ES_sector_wise_J_m1.insert(energies_ES_sector_wise_J_m1.end(),temp_Jm1);        
    }
}


M f_B_r(int r, double K, vector<vector<int>> HS_sector_1, vector<vector<int>> HS_sector_2){
    /* Calculates matrix representation of f_r^B(k) between two total momentum sectors.*/
    int dim_sector_1 = HS_sector_1.size();
    int dim_sector_2 = HS_sector_2.size();
    
    // initialize matrix to store results
    M result = M::Zero(dim_sector_1,dim_sector_2);        

    for (int l_1 = 0; l_1 < dim_sector_1; l_1++){        
        for (int l_2 = 0; l_2 < dim_sector_2; l_2++){
            vector<int> beta_1 = HS_sector_1[l_1];
            vector<int> beta_2 = HS_sector_2[l_2];
            
            result(l_1, l_2) = f_B_r_matrix_element(r, K, beta_1, beta_2); 
            }
        }
    return result;
};


double f_B_r_matrix_element(int r, double K, vector<int> beta_1, vector<int> beta_2){
    /* Calculates the matrix element of the bosonic factor f_B_r(k) (the dependence 
     * on the external momentum k is only in the total momentum sectors between 
     * which we evaluate f_B_r(k)).*/
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
        int sign = 1;
        if (abs_mi(a_m_b2) % 2 == 1){
            sign = -1;
        }
        // calculate term in sum
        result += sign * power_l_over_sqrt_Kl_mi(r, K, two_a_m_b1_m_b2) / (factorial_mi(a_m_b1) 
                * factorial_mi(a_m_b2) * factorial_mi(b1_p_b2_m_a)); 
        
        alpha = next_val(lower, upper, alpha);
    };    
    
    return sqrt(factorial_mi(beta_1)) * sqrt(factorial_mi(beta_2)) * result;
}


complex<double> K_L_summation(int k, int p_c, int r_1, int r_2, double eta, 
        double omega, eigenstates_J_pm* ES){
    /* Summation appearing in K-L-representation for G(r_1, r_2, k)*/
    complex<double> result_J_1 = 0;
    complex<double> result_J_m11 = 0;
    
    // l_1 and l_2 denote the indices of the total momentum sectors between which 
    // we calculate the matrix elements of f_B_r(k)
    for (int l_2 = min(0, -k); l_2 < 2 * p_c + 1; l_2 ++){
        int l_1 = l_2 - k;
        
        // Get eigenstates and energies of both momentum sectors
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_1_l1 = ES -> energies_ES_sector_wise_J_1[l_1];
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_1_l2 = ES -> energies_ES_sector_wise_J_1[l_2];
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_m1_l1 = ES -> energies_ES_sector_wise_J_m1[l_1];
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_m1_l2 = ES -> energies_ES_sector_wise_J_m1[l_2];
        
        int dim_sec_l1 = Eigen_J_1_l1.eigenvalues().size();
        int dim_sec_l2 = Eigen_J_1_l2.eigenvalues().size();
        
        for (int j_1 = 0; j_1 < dim_sec_l1;  j_1++){
            for (int j_2 = 0; j_2 < dim_sec_l2;  j_2++){
                double E_j1 = 0;
            }
        }
        
    }
    return 0;
}