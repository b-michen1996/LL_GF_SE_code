/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.cc to edit this template
 */

#include <iostream>
#include <cmath>


#include "multi_index_aux.h"

int momentum(int l, int p_c){
    /* Returns (integer) momentum for given index l of a mutli-index
     * and p_c */    
    int p_l =  l - p_c;
            
    if (l > p_c - 1){
        p_l = l - p_c + 1;             
    }
    
    return p_l;
}

int momentum_mi(vector<int> alpha){
    /* Returns (integer) momentum for given multi-index alpha*/    
    int p_c = int(alpha.size()/2);
    
    int p_tot = 0;     
    // calculate total momentum
    for (int l = 0; l < alpha.size(); l++ ){
        // momentum of current index
        int p_l =  momentum(l, p_c); 
        // add to total momentum
        p_tot = p_tot + p_l * alpha[l];                   
    }
    return p_tot;
}


int factorial(int l){
    /* Compute factorial*/
    int result  = 1;
    
    for (int j = 2; j < l + 1; j++){
        result = result * j;
    }
    return result;
}


int factorial_mi(vector<int> alpha){
    /* Calculate factorial of multi-index*/    
    int result = 1;
    
    for (int l = 0; l < alpha.size(); l++){
        result = result * factorial(alpha[l]);
    }
    return result;
}


int abs_mi(vector<int> alpha){
    /* Compute |alpha| for multi-index alpha. */
    int result = 0;
    
    for (int l = 0; l < alpha.size(); l++){
        result += abs(alpha[l]);
    }
    return result;
}


double energy_H0(vector<int> alpha){
    /* Energy of free Hamiltonian of a multi-index Fock state (up to prefactor).*/
    double result = 0;
    int p_c = int (alpha.size()/2 + 0.1);
    
    for (int l = 0; l < 2 * p_c; l++){
        result += abs(momentum(l, p_c)) * alpha[l];
    }
    return result;
}


double energy_H0_reg(double L, vector<int> alpha, double gamma){
    /* Energy of free Hamiltonian of a multi-index Fock state (up to prefactor).*/
    double result = 0.;
    int p_c = int (alpha.size()/2 + 0.1);
    
    double prefactor = gamma * 2 * M_PI / L;
    
    for (int l = 0; l < 2 * p_c; l++){
        result += exp(-prefactor * abs(momentum(l, p_c)))
        * abs(momentum(l, p_c)) * alpha[l];
    }
    return result;
}

vector<int> next_val(vector<int> lower, vector<int> upper, vector<int> val){
    /* return next value of multi-index val with bounds given by lower and 
     * upper. */
    int dim = lower.size();
    
    for (int l = dim - 1; l > -1; l--){
        // increment val[l] by one
        val[l] += 1;
        
        // if it exceeds the maximum value, set it to lower limit and move to next entry
        // if it does not exceed the maximum value, return it
        if (val[l] > upper[l]){
            val[l] = lower[l];
        } else{
            return val;
        }
    }
    // if all elements of the multi-index have reached the max, return empty vector
    vector<int> empty;
    return empty;
};


vector<int> next_val_Ec(vector<int> val, int E_c){
    /* return next value of multi-index with upper energy bound E_c*/
    int p_c = int(val.size()/2 +0.1);
    
    for (int l = 2 * p_c - 1; l > -1; l--){
        // increment val[l] by one
        val[l] += 1;
        
        // if energy exceeds the maximum value, set entry to zero and move to next entry
        // if it does not exceed the maximum value, return it
        if (energy_H0(val) > E_c){
            val[l] = 0;
        } else{
            return val;
        }
    }
    // if all elements of the multi-index have reached the max, return empty vector
    vector<int> empty;
    return empty;
};


double function_A(vector<int> alpha, double factor, int p_c){
    /* Function of multi-index that appears in matrix element of H1_B*/
    double result = 0;   
    
    for (int l = 0; l < 2 * p_c; l++ ){
        // momentum of current index
        int p_l =  momentum(l, p_c);        
        
        result += p_l * cos(factor * p_l) * alpha[l];        
        }
    return result;
}


double power_sqrt_l_over_l_mi(vector<int> alpha){
    /* Expression (\frac{\sqrt{l}}{l})^\alpha that appears in matrix element 
     * of H1_B*/
    double result = 1;
    int p_c = int(alpha.size() / 2);
    
    for (int l = 0; l < 2 * p_c; l++ ){
        // momentum of current index
        int p_l =  momentum(l, p_c);        
        
        result = result * pow(sqrt(abs(p_l))/p_l,  alpha[l]);        
    }
    return result;
}


double power_sqrt_l_over_l_mi_reg(double L, vector<int> alpha, double gamma){
    /* Expression (\frac{\sqrt{l}}{l})^\alpha that appears in matrix element 
     * of H1_B*/
    double result = 1;
    double prefactor = gamma * M_PI / L;
    
    int p_c = int(alpha.size() / 2);
    
    for (int l = 0; l < 2 * p_c; l++){
        // momentum of current index
        int p_l =  momentum(l, p_c);        
        
        result = result * pow(exp(-prefactor * abs(momentum(l, p_c))) 
                * sqrt(abs(p_l))/p_l ,  alpha[l]);        
    }
    return result;
}


double power_l_over_sqrt_Kl_mi(int r, double K, vector<int> alpha){
    /* Expression that appears in matrix element 
     * of f_r^B(k)*/
    double result = 1;
    int p_c = int(alpha.size() / 2 + 0.1);
    
    for (int l = 0; l < 2 * p_c; l++ ){
        // momentum of current index
        int p_l =  momentum(l, p_c);        
        int sign = 1.;
        if (p_l < 0){
            sign = -1.;
        }
        if (alpha[l] > 0){
            result = result * pow((1 + sign *  r *  K) / (2 * sqrt(K * abs(p_l))),  alpha[l]);
        };
        
    }
    return result;
}


double power_l_over_sqrt_Kl_mi_reg(int r, double K, double L, double gamma, vector<int> alpha){
    /* Expression that appears in matrix element 
     * of f_r^B(k)*/
    double result = 1;
    double prefactor = gamma * M_PI / L;
    
    int p_c = int(alpha.size() / 2 + 0.1);
    
    for (int l = 0; l < 2 * p_c; l++ ){
        // momentum of current index
        int p_l =  momentum(l, p_c);        
        int sign = 1;
        if (p_l < 0){
            sign = -1;
        }
        if (alpha[l] > 0){
        result = result * pow(exp(-prefactor * abs(momentum(l, p_c))) 
                * (1 + sign *  r *  K) / (2 * sqrt(K * abs(p_l))),  alpha[l]);
        };
        
    }
    return result;
}