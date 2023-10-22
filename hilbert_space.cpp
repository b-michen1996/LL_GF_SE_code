/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */

#include "hilbert_space.h"
#include "multi_index_aux.h"

#include <cmath>
#include <iostream>


using namespace std;

HS::HS(int in_pc):p_c(in_pc){  
    /* Constructor for HS struct containing all sectors of total Hilbert space*/
    HS_tot.resize(2 * p_c + 1);
    
    vector<int> state;
    fill_state_list(HS_tot, state, 0, p_c);           
    
};

void fill_state_list(vector<vector<vector<int>>>&  list, vector<int> state_curr, int E_curr, int p_c){
    /* Recursive function to loop over all allowed states in the truncated Hilbert
     * space and assign them to their respective total momentum sector */
    
    // check if we reached the full length, if yes calculate total momentum and assign 
    // to sector of Hilbert space
    int state_length = state_curr.size();

    if (state_length == 2 * p_c ){
        int p_tot = 0;     
        // calculate total momentum
        for (int l = 0; l < 2 * p_c; l++ ){
            // momentum of current index
            int p_l =  momentum(l, p_c); 
            
            p_tot = p_tot + p_l * state_curr[l];                   
        }
        
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