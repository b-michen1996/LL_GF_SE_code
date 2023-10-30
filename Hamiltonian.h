/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   H_1_sector.h
 * Author: bm
 *
 * Created on October 18, 2023, 4:57 PM
 */

#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>
#include <array>

#include "hilbert_space.h"

using namespace std;
typedef Eigen::MatrixXcd M;


M H_Luttinger_J(double v_F, double K, double U_m, double L, double a, double alpha, vector<vector<int>> HS_sector, int J);

double H1_B_matrix_element(double K, vector<int> beta_1, vector<int> beta_2, double L, double a);

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
