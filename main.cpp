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

#include "hilbert_space.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {   
    // some parameters
    double v_F;
    double L = 100;
    double a = 0.1;        
    double E_c = 10;
    
    int run_Nr = 1;
    
    int p_c = 40;
    int index = 3;
    
    
    HS hd_test(p_c);
    
    
    cout << hd_test.HS_tot.size() << "\n";
    
    /*
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
    }
    */
    
    

    
    

    return 0;
}

