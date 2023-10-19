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
    
    HS_tot.resize(2 * p_c);
    
    vector<int> state;
    fill_state_list(HS_tot, state, 0, p_c);        
            
         
    
};

void fill_state_list(vector<vector<vector<int>>>&  list, vector<int> state_curr, int E_curr, int p_curr, int p_c){
    /* Recursive function to loop over all allowed states in the truncated Hilbert
     * space and assign them to their respective total momentum sector */
    
    // check if we reached the full length, calculate total momentum and assign 
    // to sector of Hilbert space
    int state_length = state_curr.size();
    if (state_length == 2 * p_c ){
        int p_tot = 0;
        
        // calculate total momentum
        for (int l = -p_c; l < p_c + 1; l++ ){
            // skip l = 0
            if (l == 0){
            continue;
            }        
            p_tot = l * state_curr[l + p_c];        
        }
        
        // index for total momentum sector
        int l_p_tot = 0;
                
        if (p_tot < 0){ 
            l_p_tot = p_tot + p_c;
        }else{
            l_p_tot = p_tot + p_c - 1;
        }
                    
        // insert vector into list of states for respective total momentum sector
        list[l_p_tot].insert(list[l_p_tot].end(), state_curr);
    }else {
        for (int l = 0; l++; l < p_c - E_curr){
            int E = E_curr + l * state_length;
        }
        
    }

    return;
}