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
    int p_c;
    vector<vector<vector<int>>>  HS_tot;
    vector<vector<int>> index_occupation_number;
        
    //constructor method
    HS(int in_pc);

    // return vector of Fock states (each in the form of a vector of occupation 
    // numbers) in a given total momentum sector up to cut-off energy E_c
    vector<vector<int>> sector(int p);
};


struct HS_Ec_pc{
    /* Struct that will create list of all Fock states belonging to the individual total 
     * momentum sectors of the Hilbert space under constrained by the cut-off energy
     * E_c (here taken as integer ignoring all prefactors). The momentum cut-off 
     * p_c is independent. */
    int E_c;
    int p_c;
   
    vector<vector<vector<int>>>  HS_tot;
            
    //constructor method
    HS_Ec_pc(int in_Ec, int in_pc);
};


void fill_state_list(vector<vector<vector<int>>>&  list, vector<int> state_curr, int E_curr,  int p_c);

