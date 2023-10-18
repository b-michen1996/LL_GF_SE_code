/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   hilbert_space.h
 * Author: bm
 *
 * Created on October 18, 2023, 9:41 PM
 */

#pragma once
#include <vector>

using namespace std;

struct HS{
    /* Struct that will create list of all Fock states belonging to the individual total 
     * momentum sectors of the Hilbert space under consideration of the cut-off energy
     * E_C = v_F |P_c|. Everything is entirely determined by the positive integer 
     * p_c = |P_c| L / (2pi) */
    int pc;
    vec<vec<vec<int>>>  HS_tot;
    
    int Ny;
    double Z;
    double deltaNL;
    double c;
    double delta;  
    /*
    Eigen::SparseMatrix<std::complex<double>> H1_stat;
    Eigen::SparseMatrix<std::complex<double>> H2_stat;
    Eigen::SparseMatrix<std::complex<double>> H3_stat;
    Eigen::SparseMatrix<std::complex<double>> H4_stat;
     */
    
HS(int in_pc);

//void operator()(const state_type& x, state_type& dxdt, double t);
};


