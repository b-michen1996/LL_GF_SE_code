/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */

#include "hilbert_space.h"
#include <cmath>
#include <iostream>

using namespace std;

vector<int> occupation_number(int index, int p_c){
    /* Returns vector of occupation numbers for state with given 
     * index on truncated Hilbert space */
    
    // total number of states in truncated HS
    int n_tot = int(pow(tgamma(p_c + 1), 2));
        
    vector<int> result;
        
    // fill vector of occupation numbers 
    for (int l = p_c; -p_c <= l; l-- ){
        // skip l = 0
        if (l == 0){
            continue;
        }
        
        // number of states up to index of given mode
        int A = 0;
        if(l > 0){
            A =  tgamma(p_c - l + 1);
        }else{
            A =  tgamma(p_c + 1) * tgamma(abs(l) + 1);
        }  
                        
        // maximum occupation number 
        int n_l_max = int(abs(p_c / l));  
            
        // occupation number for given mode             
        int n_l = int (index / A) % n_l_max;
            
        result.insert(result.end(), n_l);  
    }    
    return result; 
    
}

HS::HS(int in_pc):p_c(in_pc){  
    /* Constructor for HS struct containing all sectors of total Hilbert space*/
    int n_tot = int(pow(tgamma(p_c + 1), 2));
    
    
    
    for (int j = 0; j++; j <= n_tot){
        // fill vector of occupation number / index mapping
        for (int l = p_c; l--; -p_c <= l ){
            // occupation number for given mode 
            int n_l = j % int((tgamma(p_c - l + 1)));
        }
    }
    
    for (int j = p_c; j--; -p_c <= j ){
            // maximum occupation number 
            int n_j_max = int(abs(p_c / j));
            
            
            
        }
    
    
};
