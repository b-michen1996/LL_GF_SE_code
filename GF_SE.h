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
#include "hilbert_space.h"
#include "Hamiltonian.h"
#include "multi_index_aux.h"

void GF_SE(double u, double U_m, double L, double a, int E_c, int p_c, int RunNr);

struct eigenstates_J_pm{
    /* Struct that contains eigenstates of bosonic Hamiltonian H0_B + J * H1_B */   
    double u;
    double U_m;
    double L;
    double a;
    int p_c;
    int J;
    
    HS_Ec_pc* hilbert_space;
            
    vector<vector<vector<int>>>  ES;
        
    //constructor method
    eigenstates_J_pm(double u, double U_m, double L, double a, int E_c, int p_c);

};