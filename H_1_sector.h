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

using namespace std;
typedef Eigen::MatrixXcd M;

M H1_B(double L, double a, vector<vector<int>> HS_sector);

double H1_B_matrix_element(vector<int> beta_1, vector<int> beta_2, double L, double a);

