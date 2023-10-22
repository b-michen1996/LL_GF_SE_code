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

#include "hilbert_space.h"
#include "H_1_sector.h"
#include "multi_index_aux.h"

using namespace std;
typedef Eigen::MatrixXcd M;
/*
 * 
 */
int main(int argc, char** argv) {   
    // some parameters
    double v_F;
    double L = 100;
    double a = 1;        
    double E_c = 10;
    
    int run_Nr = 1;
    
    int p_c = 15;
    int index = 3;
    
    vector<int> beta_1 = {1,1,2,3,0,1};
    vector<int> beta_2 = {1,2,2,3,0,5};
    
    /*
    double ma = H1_B_matrix_element(beta_1, beta_2, L, a);
    
    cout << ma;
    HS hd_test(p_c);
     */ 
    
    HS hd_test(p_c);
        
    for (int l = 0; l < 2 * p_c + 1; l++){ 
        vector<vector<int>> m_sec_curr = hd_test.HS_tot[l];
        M H1_block = H1_B(L, a, m_sec_curr);           
    }
    
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
    
    
    

    
    

    return 0;
}

