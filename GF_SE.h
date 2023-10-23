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

void GF_SE(double L, double a, int p_c);

struct eigenstates_J_pm{
    /* Struct that contains eigenstates of bosonic Hamiltonian H0_B + J * H1_B */
   
    double L;
    double a;
    int p_c;
    int J;
    
    HS* hilbert_space;
            
    vector<vector<vector<int>>>  eigenstates;
        
    //constructor method
    eigenstates_J_pm(double int in_pc);

};