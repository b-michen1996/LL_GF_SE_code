/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */

#include "H_1_sector.h"
#include <cmath>
#include <array>
#include <iostream>
#include <omp.h>



M H1_B(double L, double a, vector<vector<int>> HS_sector){
    /* Calculates matrix representation of bosonic Hamiltonian H1_B on the subspace given by 
     * the Fock states / vectors of occupation numbers in HS_sector. */
    int dim_sector = HS_sector.size();
    
    // initialize matrix to store results
    M result = M::Zero(dim_sector,dim_sector);
    
    
    // initialize empty sparse matrix H1_B
    
    for (int l_1 = 0; l_1 < dim_sector; l_1++){
        for (int l_2 = 0; l_2 < dim_sector; l_2++){
            }
        }
    
    return result;
};

double matrix_element(vector<int> beta_1, vector<int> beta_2, double L, double a){
    /* Calculate matrix element of H1_B between two Fock states beta_1, beta_2*/
    int dim = beta_1.size();
    
    // We need to loop over multi-index alpha, this gives array of lower and 
    // upper bounds for each component.
    vector<int> lower(dim);
    vector<int> upper(dim);
    
    for (int l= 0; l < dim; l++){  
        lower[l] = max(beta_1[l], beta_2[l]);
        upper[l] = beta_1[l] + beta_2[l];
    }
    
    // loop over all allowed values of alpha using auxilliary function
    
    return 0;
}

vector<int> next_val(vector<int> lower, vector<int> upper, vector<int> val){
    /* return next value of multi-index val with bounds given by lower and upper. */
};