/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   Export_MA.h
 * Author: bm
 *
 * Created on November 3, 2023, 8:50 AM
 */

#pragma once

#include "GF_SE.h"

void export_f_B(double v_F, double K, double U_m, double L, double a, double alpha, 
        int E_c, int p_c, int runNr, int threads);

void write_f_B_to_file(M f_B_R_curr_ES_base, string filename);