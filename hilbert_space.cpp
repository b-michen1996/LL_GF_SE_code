/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */

#include "hilbert_space.h"
#include "multi_index_aux.h"

#include <cmath>
#include <iostream>
#include <array>


using namespace std;

HS::HS(int in_pc):p_c(in_pc){  
    /* Constructor for HS struct containing all sectors of total Hilbert space*/
    HS_tot.resize(2 * p_c + 1);     
    
    vector<int> state;
    fill_state_list(HS_tot, state, 0, p_c);           
    
};


HS_Ec_pc::HS_Ec_pc(int in_Ec, int in_pc):E_c(in_Ec), p_c(in_pc){  
    /* Constructor for HS_Ec_pc struct containing all sectors of total Hilbert 
     * space that is truncated up to single-particle momentum p_c and energy E_c*/
    HS_tot.resize(2 * E_c + 1);     
    
    // starting state of zeros
    vector<int> state(2 * p_c, 0);
    
    // loop over all states up to energy cut-off E_c
    while (state.size() > 0){
        // calculate total momentum
        int p_tot = momentum_mi(state);
        
        // index for total momentum sector (total momenta run from -E_c 
        // to E_c and not -p_c to p_c!!
        int l_p_tot = p_tot + E_c;
                                    
        // insert vector into list of states for respective total momentum sector
        HS_tot[l_p_tot].insert(HS_tot[l_p_tot].end(), state);
        
        state = next_val_Ec(state, E_c);
    }       
};


void fill_state_list(vector<vector<vector<int>>>&  list, vector<int> state_curr, int E_curr, int p_c){
    /* Recursive function to loop over all allowed states in the truncated Hilbert
     * space and assign them to their respective total momentum sector */
    
    // check if we reached the full length, if yes calculate total momentum and assign 
    // to sector of Hilbert space
    int state_length = state_curr.size();
    
    if (state_length == 2 * p_c ){
        int p_tot = momentum_mi(state_curr);
        
        // index for total momentum sector
        int l_p_tot = p_tot + p_c;
                                    
        // insert vector into list of states for respective total momentum sector
        list[l_p_tot].insert(list[l_p_tot].end(), state_curr);
        
        
    }else { 
        int pl_curr = momentum(state_length, p_c);
        
        for (int n = 0; n * abs(pl_curr)  < p_c - E_curr +1; n++){
             
            int E = E_curr + n * abs(pl_curr);            
         
            vector<int> new_state = state_curr;           
            new_state.insert(new_state.end(), n);            
            fill_state_list(list, new_state, E, p_c);            
        }        
    }

    return;
}


