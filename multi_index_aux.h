/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   multi_index_aux.h
 * Author: bm
 *
 * Created on October 22, 2023, 6:05 PM
 */

#pragma once

#include <vector>


using namespace std;

int momentum(int l, int p_c);

int momentum_mi(vector<int> alpha);

int factorial(int l);

int factorial_mi(vector<int> alpha);

int abs_mi(vector<int> alpha);

double energy_H0(vector<int> alpha);

vector<int> next_val(vector<int> lower, vector<int> upper, vector<int> val);

vector<int> next_val_Ec(vector<int> val, int E_c);

double function_A(vector<int> alpha, double factor, int p_c);

double power_sqrt_l_over_l_mi(vector<int> alpha);