//
//  support.cpp
//  travelling_wave
//
//  Created by Mirco Meazzo on 11/09/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#include "support.hpp"
#include <iostream>


void Grid::initialize_grid(){
    
    double step = length / nx;
#pragma omp parallel for
    for (int i=0; i<nx; i++) {
        position[i]=step*i;
    }
}



void Grid::show_grid(){
    
    using namespace std;
    for (int i=0; i<nx; i++) {
        cout << "[" << i << "]= " << this->position[i] << endl;
    }
}



double* Grid::compute_deltax(){
    
    static double* deltax= new double(nx);
    
    for (int i=0; i<nx; i++) {
        deltax[i] = this->position[i+1] - this->position[i];
    }
    
    return deltax;
}



double* f::solve_pde(double c) {
    
//    build rhs
    double* rhs = new double(nx);
    
//    centered FD
    for (int i=1; i<nx-1; i++) {
        rhs[i]= f_t[i-1] + f_t[i+1] - 2*f_t[i];
    }
//    BCs
    rhs[0] = f_t[0]+f_t[2]- 2*f_t[1];
    rhs[nx]= f_t[nx]+f_t[nx-2] -2*f_t[nx-1];
    
    for (int i=0; i<nx; i++) {
        f_tn[i]= c*c*rhs[i] + 2*f_t[i]- f_tm[i];
    }
    return f_tn;
}
