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

using namespace std;
using namespace std::chrono;

typedef Eigen::MatrixXcd M;
/*
 * 
 */
int main(int argc, char** argv) {   
    // some parameters
    double v_F;
    double L = 100;
    double a = 1;        
    double u = 1.;
    double U_m = 0.1;
    int E_c = 15;
    
    int run_Nr = 1;
    
    int p_c = 14;
    int index = 3;
    
    /*
    auto t1 = std::chrono::system_clock::now();
    eigenstates_J_pm diagonalization_test(u, U_m, L, a, E_c, p_c);
    auto t2 = high_resolution_clock::now();  
    auto duration = duration_cast<milliseconds>(t2 - t1);
    
    cout << "Duration of full diagonalization " << duration.count() <<"ms \n";
    */
    
    auto t1 = std::chrono::system_clock::now();
    HS HS_test(p_c);
    auto t2 = high_resolution_clock::now();  
    auto duration = duration_cast<milliseconds>(t2 - t1);
    cout << "Duration of HS " << duration.count() <<"ms \n";
    
    t1 = std::chrono::system_clock::now();
    HS_Ec_pc HS_Ec_pc_test(E_c, p_c);
    t2 = high_resolution_clock::now();  
    duration = duration_cast<milliseconds>(t2 - t1);
    cout << "Duration of HS_Ec_pc " << duration.count() <<"ms \n";
    
    for (int l = 0; l < 2 * p_c + 1; l++){   
        int counter = 0;
        vector<vector<int>> m_sec_curr_HS = HS_test.HS_tot[l];
        vector<vector<int>> m_sec_curr_HS_Ec_pc = HS_Ec_pc_test.HS_tot[l];
        
        cout << "l = " << l << ", size of sector " << m_sec_curr_HS.size() << " vs"
                << m_sec_curr_HS_Ec_pc.size() << "\n";
        cout << "--------------------------------------\n";
        
        for (int j = 0; j < m_sec_curr_HS_Ec_pc.size(); j++){
            vector<int> ES_1 = m_sec_curr_HS[j];
            vector<int> ES_2 = m_sec_curr_HS_Ec_pc[j];
            
            for (int k = 0; k < 2 * p_c; k++){
                if (abs(ES_1[k] - ES_2[k]) > 0){
                    cout << "Attention";
                }               
            }
        }
        cout << "--------------------------------------\n";
        
        /*
        for (vector<int> oc_n : m_sec_curr_HS){
            //counter += 1;
            //cout << "state number " << counter << " in sector " << l << "\n";
            for (int n : oc_n){                   
                cout << n << ", ";
            }
            cout << "\n";            
        }
        
        for (vector<int> oc_n : m_sec_curr_HS_Ec_pc){
            //counter += 1;
            //cout << "state number " << counter << " in sector " << l << "\n";
            for (int n : oc_n){                   
                cout << n << ", ";
            }
            cout << "\n";            
        }
        cout << "--------------------------------------\n";*/
    }
    
    return 0;
}
    /* 
    vector<int> beta_1 = {0,2,0,0,0,0};
    
    vector<int> beta_2 = {0,0,0,0,2,0};
    
    beta_1 = {3,0,0,0,5,5};
    
    beta_2 = {6,0,0,2,2,6};
    

    
    
    HS hd_test(p_c);
    
    vector<int> beta_1 = hd_test.HS_tot[0][0];
    vector<int> beta_2 = hd_test.HS_tot[0][int(p_c/4)];
   
    auto t1 = std::chrono::system_clock::now();
    
    double ma = H1_B_matrix_element(beta_1, beta_2, L, a);
    
    auto t2 = high_resolution_clock::now();  
    auto duration = duration_cast<nanoseconds>(t2 - t1);
    cout << "Matrix element " << ma << ", duration " << duration.count() <<"ns \n";
     */
    /*
    cout << ma;
    HS hd_test(p_c);
     */ 
    
    /*
    HS hd_test(p_c);        
    
    for (int l = 0; l < 2 * p_c + 1; l++){          
        vector<vector<int>> m_sec_curr = hd_test.HS_tot[l];
        
        auto t1 = std::chrono::system_clock::now();
        
        M H1_block = H1_B(L, a, m_sec_curr);           
        auto t2 = high_resolution_clock::now();
        
        auto duration = duration_cast<milliseconds>(t2 - t1);
        
        cout << "Current momentum sector " << l - p_c << ", size of sector " << m_sec_curr.size() <<"\n";
        cout << ", duration " << duration.count() <<"ms \n";
    }
    */
    /*
    HS hd_test(p_c);
    
    
    cout << hd_test.HS_tot.size() << "\n";
    
    
    for (int l = 0; l < 2 * p_c + 1; l++){   
        int counter = 0;
        vector<vector<int>> m_sec_curr = hd_test.HS_tot[l];
        
        cout << "l = " << l << ", size of sector " << m_sec_curr.size() << "\n";
        
        for (vector<int> oc_n : m_sec_curr){
            counter += 1;
            cout << "state number " << counter << " in sector " << l << "\n";
            for (int n : oc_n){                   
                cout << n << ", ";
            }
            cout << "\n";            
        }
    }*/
    
    
    

    
    



