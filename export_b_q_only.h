/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   export_bq_only.h
 * Author: bm
 *
 * Created on November 10, 2023, 4:19 PM
 */

#pragma once

#include "GF_SE.h"

void export_b_q_only(double v_F, double K, double U_m, double L, double a, double gamma, 
        int E_c, int p_c, int runNr, int threads);

M b_q_Fock(int l, vector<vector<int>> HS_sector_1, vector<vector<int>> HS_sector_2);



