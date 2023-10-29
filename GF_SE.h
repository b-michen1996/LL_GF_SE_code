/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   GF_SE.h
 * Author: bm
 *
 * Created on October 23, 2023, 12:00 PM
 */

#pragma once

#include <vector>
#include <eigen3/Eigen/Eigenvalues>
#include <cmath>


#include "hilbert_space.h"
#include "Hamiltonian.h"

void GF_SE_J_explicit(double v_F, double K, double U_m, double L, double a, double alpha, double eta, 
        int E_c, int p_c, double omega, double beta, int runNr, int threads);

struct eigenstates_J_pm{
    /* Struct that contains eigenstates of bosonic Hamiltonian H0_B + J * H1_B
     * for both values of J = 1, -1.*/   
    double v_F;
    double K;
    double U_m;
    double L;
    double a;
    double alpha;
    int E_c;
    int p_c;
    int J;
    int threads;
    
    // Declare member of struct containing truncated Hilbert space sectors
    HS_Ec_pc HS_truncated{E_c, p_c};
    
    // Vectors containing Eigensolvers filled with eigenstates and energies of
    // each sector
    vector<Eigen::SelfAdjointEigenSolver<M>>  energies_ES_sector_wise_J_1;
    vector<Eigen::SelfAdjointEigenSolver<M>>  energies_ES_sector_wise_J_m1;
            
    //constructor method
    eigenstates_J_pm(double in_v_F, double in_K, double in_U_m, double in_L, double in_a, 
    double in_alpha, int in_E_c, int in_p_c, int in_threads);
};


M f_B_r(int r, double K, vector<vector<int>> HS_sector_1, vector<vector<int>> HS_sector_2);

double f_B_r_matrix_element(int r, double K, vector<int> beta_1, vector<int> beta_2);


complex<double> K_L_summation(int k, int p_c, int r_1, int r_2, double eta, 
        double omega, eigenstates_J_pm* ES);
 