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
typedef Eigen::MatrixXd M;

M H1_B(double L, double a, int P_tot, int pc, vector<vector<int>> HS_sector);

double matrix_element(vector<int> beta_1, vector<int> beta_2, double L, double a);

vector<int> next_val(vector<int> lower, vector<int> upper, vector<int> lower val);
