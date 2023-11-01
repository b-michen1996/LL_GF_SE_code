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
    ofstream parameters(output_folder + "/parameters.csv");    
    
    parameters << "v_F " << v_F << "\n" << "K " << K << "\n" << "alpha " << alpha << "\n" 
            << "eta " << eta << "\n" << "U_m " << U_m << "\n" << "L " << L << "\n" 
            << "a " << a << "\n" << "E_c " << E_c << "\n" << "p_c " << p_c << "\n"  
            << "omega " << omega << "\n" << "beta " << beta << "\n" 
            << "runNr " << runNr << "\n" << "threads " << threads;
    parameters.close();
    
    // Calculate eigenstates
    eigenstates_J_pm ES(v_F, K, U_m, L, a, alpha, E_c, p_c, threads);
        
    // Calculate partition function
    double Z = 0;
    
    for (int l = 0; l < 2 * E_c + 1; l++){
        double Z_sector_J_1 = 0;
        double Z_sector_J_m1 = 0;
        
        // get dimension of current momentum sector
        int dim_sector = ES.HS_truncated.HS_tot[l].size();
        // get current sector
        Eigen::SelfAdjointEigenSolver<M> solver_curr_J_1 = ES.energies_ES_sector_wise_J_1[l]; 
        Eigen::SelfAdjointEigenSolver<M> solver_curr_J_m1 = ES.energies_ES_sector_wise_J_m1[l]; 
        for (int j = 0; j < dim_sector; j++){
            double E_J_1_j = solver_curr_J_1.eigenvalues().coeffRef(j);
            double E_J_m1_j = solver_curr_J_m1.eigenvalues().coeffRef(j);
            // cout << E_J_1_j << ", " << E_J_m1_j << "\n";
            Z_sector_J_1 += exp(- beta * E_J_1_j);
            Z_sector_J_m1 += exp(- beta * E_J_m1_j);
        }
        /*
        cout << "Z_sector_J_1 = " << Z_sector_J_1 << ", Z_sector_J_m1 = " 
                << Z_sector_J_m1 << "\n";
        */
        Z += Z_sector_J_1 + Z_sector_J_m1;        
    }
    
    // calculate elements of retarded GF for each k from -p_c to p_c and store 
    // in vector
    vector<complex<double>> G_RR(2 * p_c + 1, 0.);
    vector<complex<double>> G_RL(2 * p_c + 1, 0.);
    vector<complex<double>> G_LR(2 * p_c + 1, 0.);
    vector<complex<double>> G_LL(2 * p_c + 1, 0.);
    
    double prefactor = pow(2 * M_PI * alpha / L, pow(1 - K, 2) / (2 * K)) / (Z * pow(L, 2));
    
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (int l = 0; l < 2 * p_c + 1; l++ ){
        int k = l - p_c;
        int r_1 = 1;
        int r_2 = -1;
        
        vector<complex<double>> G_R_R = K_L_summation(k, E_c, 1, 
                1, K, eta, omega, beta, &ES);
        G_RR[l] = prefactor * (G_R_R[0] + G_R_R[1]);
        /*
        cout << "k = " << k << ", G_R_R = " << prefactor * G_R_R[0] << ", " << prefactor * G_R_R[1] 
                << prefactor * G_R_R[2] << ", " << prefactor * G_R_R[3] << "\n";
        */
        
        vector<complex<double>> G_R_L = K_L_summation(k, E_c, 1, 
                -1, K, eta, omega, beta, &ES);
        G_RL[l] = prefactor * (G_R_L[0] + G_R_L[1]);
        /*
        cout << "k = " << k << ", G_R_L = " << G_R_L[0] << ", " << G_R_L[1] 
                << G_R_L[2] << ", " << G_R_L[3] << "\n";
        */
        vector<complex<double>> G_L_R = K_L_summation(k, E_c, -1, 
                1, K, eta, omega, beta, &ES);
        G_LR[l] = prefactor * (G_L_R[0] + G_L_R[1]);
        /*
        cout << "k = " << k << ", G_L_R = " << G_L_R[0] << ", " << G_L_R[1] 
                << G_L_R[2] << ", " << G_L_R[3] << "\n";
        */
        vector<complex<double>> G_L_L = K_L_summation(k, E_c, -1, 
                -1, K, eta, omega, beta, &ES);
        G_LL[l] = prefactor * (G_L_L[0] + G_L_L[1]);
        /*
        cout << "k = " << k << ", G_L_L = " << G_L_L[0] << ", " << G_L_L[1] 
                << G_L_L[2] << ", " << G_L_L[3] << "\n";
        */
        cout << "Calculation done for k = " << k << "\n";
    } 
    
    // save result to file
    ofstream G_RR_file(output_folder + "/G_RR.csv");
    ofstream G_RL_file(output_folder + "/G_RL.csv");
    ofstream G_LR_file(output_folder + "/G_LR.csv");
    ofstream G_LL_file(output_folder + "/G_LL.csv");
    
    for (int l = 0; l < 2 * p_c + 1; l++ ){
        if (G_RR[l].imag() > 0){
            G_RR_file << G_RR[l].real() << "+" << G_RR[l].imag() << "j \n";
        }
        else{
            G_RR_file << G_RR[l].real() << "-" << abs(G_RR[l].imag()) << "j \n";
        } 
        
        if (G_RL[l].imag() > 0){
            G_RL_file << G_RL[l].real() << "+" << G_RL[l].imag() << "j \n";
        }
        else{
            G_RL_file << G_RL[l].real() << "-" << abs(G_RL[l].imag()) << "j \n";
        } 
        
        if (G_LR[l].imag() > 0){
            G_LR_file << G_LR[l].real() << "+" << G_LR[l].imag() << "j \n";
        }
        else{
            G_LR_file << G_LR[l].real() << "-" << abs(G_LR[l].imag()) << "j \n";
        } 
        
        if (G_LL[l].imag() > 0){
            G_LL_file << G_LL[l].real() << "+" << G_LL[l].imag() << "j \n";
        }
        else{
            G_LL_file << G_LL[l].real() << "-" << abs(G_LL[l].imag()) << "j \n";
        } 
    }
    G_RR_file.close();
    G_RL_file.close();
    G_LR_file.close();
    G_LL_file.close();
    return;
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


vector<complex<double>> K_L_summation(int k, int E_c, int r_1, int r_2, double K, double eta, 
        double omega, double beta, eigenstates_J_pm* ES){
    /* Summation appearing in K-L-representation for G(r_1, r_2, k)*/
    complex<double> result_J_1 = 0.;
    complex<double> result_J_m1 = 0.;
    complex<double> result_only_1 = 0.;
    complex<double> result_only_m1 = 0.;

    // l_1 and l_2 denote the indices of the total momentum sectors between which 
    // we calculate the matrix elements of f_B_r(k)
    for (int l_2 = max(0, k); l_2 < 2 * E_c + 1; l_2 ++){
        int l_1 = l_2 - k;
        if (l_1 > 2 * E_c){
            break;
        }
        
        // Get eigenstates and energies of both momentum and J sectors 
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_1_l1 = ES -> energies_ES_sector_wise_J_1[l_1];
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_1_l2 = ES -> energies_ES_sector_wise_J_1[l_2];
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_m1_l1 = ES -> energies_ES_sector_wise_J_m1[l_1];
        Eigen::SelfAdjointEigenSolver<M> Eigen_J_m1_l2 = ES -> energies_ES_sector_wise_J_m1[l_2];
        
        // get matrix rep of f_B_r(k)
        vector<vector<int>> HS_sector_l1 = ES -> HS_truncated.HS_tot[l_1];
        vector<vector<int>> HS_sector_l2 = ES -> HS_truncated.HS_tot[l_2];        
                
        M f_r1 = f_B_r(r_1, K, HS_sector_l1, HS_sector_l2);
        M f_r2 = f_B_r(r_2, K, HS_sector_l1, HS_sector_l2);
                        
        // carry out summation    
        int dim_sec_l1 = Eigen_J_1_l1.eigenvalues().size();
        int dim_sec_l2 = Eigen_J_1_l2.eigenvalues().size();
        

        for (int j_1 = 0; j_1 < dim_sec_l1;  j_1++){
            for (int j_2 = 0; j_2 < dim_sec_l2;  j_2++){  
                // translate into const to avoid conflict with eigen. This is 
                // not nice!
                const int j_1_c = j_1;
                const int j_2_c = j_2;
                
                // get eigenenergies
                double E_J_1_j1 = Eigen_J_1_l1.eigenvalues().coeffRef(j_1_c);
                double E_J_1_j2 = Eigen_J_1_l2.eigenvalues().coeffRef(j_2_c); 
                double E_J_m1_j1 = Eigen_J_m1_l1.eigenvalues().coeffRef(j_1_c);
                double E_J_m1_j2 = Eigen_J_m1_l2.eigenvalues().coeffRef(j_2_c); 
                
                // get eigenstates
                Eigen::VectorXcd EV_J_1_j1 = Eigen_J_1_l1.eigenvectors().col(j_1_c);
                Eigen::VectorXcd EV_J_1_j2 = Eigen_J_1_l2.eigenvectors().col(j_2_c);
                Eigen::VectorXcd EV_J_m1_j1 = Eigen_J_m1_l1.eigenvectors().col(j_1_c);
                Eigen::VectorXcd EV_J_m1_j2 = Eigen_J_m1_l2.eigenvectors().col(j_2_c);                
                
                result_J_1 += EV_J_1_j1.dot(f_r1 * EV_J_m1_j2) * EV_J_1_j1.dot(f_r2 * EV_J_m1_j2) 
                        * (exp(- beta * E_J_1_j1) + exp(- beta * E_J_m1_j2)) / (omega + E_J_1_j1
                        - E_J_m1_j2 + 1i * eta);
                
                result_J_m1 += EV_J_m1_j1.dot(f_r1 * EV_J_1_j2) * EV_J_m1_j1.dot(f_r2 * EV_J_1_j2) 
                        * (exp(- beta * E_J_m1_j1) + exp(- beta * E_J_1_j2)) / (omega + E_J_m1_j1
                        - E_J_1_j2 + 1i * eta);
                
                result_only_1 += EV_J_1_j1.dot(f_r1 * EV_J_1_j2) * EV_J_1_j1.dot(f_r2 * EV_J_1_j2) 
                        * (exp(- beta * E_J_1_j1) + exp(- beta * E_J_1_j2)) / (omega + E_J_1_j1
                        - E_J_1_j2 + 1i * eta);
                
                result_only_m1 += EV_J_m1_j1.dot(f_r1 * EV_J_m1_j2) * EV_J_m1_j1.dot(f_r2 * EV_J_m1_j2) 
                        * (exp(- beta * E_J_m1_j1) + exp(- beta * E_J_m1_j2)) / (omega + E_J_m1_j1
                        - E_J_m1_j2 + 1i * eta);
            }
        }
        
    }
    
    vector<complex<double>> result = {result_J_1, result_J_m1, result_only_1, result_only_m1};
    
    
    return result;
}