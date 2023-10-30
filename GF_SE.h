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

#include "Hamiltonian.h"

void GF_SE_J_explicit(double v_F, double K, double U_m, double L, double a, double alpha, double eta, 
        int E_c, int p_c, double omega, double beta, int runNr, int threads);


M f_B_r(int r, double K, vector<vector<int>> HS_sector_1, vector<vector<int>> HS_sector_2);

double f_B_r_matrix_element(int r, double K, vector<int> beta_1, vector<int> beta_2);


vector<complex<double>> K_L_summation(int k, int E_c, int r_1, int r_2, double K, double eta, 
        double omega, double beta, eigenstates_J_pm* ES);
 