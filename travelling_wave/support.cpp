//
//  support.cpp
//  travelling_wave
//
//  Created by Mirco Meazzo on 11/09/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#include "support.hpp"
#include <iostream>
#include <fstream>


void Grid::initialize_grid(){
    
    double step = length / nx;
#pragma omp parallel for schedule(dynamic)
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



double Grid::compute_deltax(){
    
    double deltax;
    deltax = this->position[1] - this->position[0];
    
    return deltax;
}



void f::init_f(int nx, double arr_t[]) {
    
    for (int i=0; i<nx; i++) {
        f_t[i]= arr_t[i];
    }
}



void f::solve_pde(double c, double t_limit) {
    
    std::ofstream outfile;
    outfile.open("results.txt", std::ios::trunc);
    if (outfile.is_open()) {
        for (int i=0; i<nx; i++) {
            outfile << f_t[i] << " ";
            if (i==nx-1) {
                outfile << std::endl;
            }
        }
        
        double t=0, dt;
        double* rhs = new double(nx);
        double dx;
        dx = f_grid.compute_deltax();
        
        dt = dx/c;
        double C = c*dt/dx;
        std::cout << "Courant number of the simulation set to " << C << std::endl;
        
        
        while (t<t_limit) {
            //    build rhs
#pragma omp parallel for schedule(dynamic)
            for (int i=0; i<nx; i++) {
                rhs[i]= f_t[i-1] + f_t[i+1] - 2*f_t[i];
            }
            if (t==0) {
            //        Solve f at 0+dt
#pragma omp parallel for schedule(dynamic)
                for (int i=0; i<nx; i++) {
                    f_tn[i]= -0.5*C*C*rhs[i] + f_t[i];
                }
            }
            else {
            //        Solve f at t+dt
#pragma omp parallel for schedule(dynamic)
                for (int i=0; i<nx; i++) {
                    f_tn[i]= C*C*rhs[i] + 2*f_t[i]- f_tm[i];
                }
            }
            f_tn[0]=f_tn[nx-1]=0;
    
    //        Update f in time
            t=t+dt;
//            Write solution to disk and move it to f_t
            for (int i=0; i<nx; i++) {
                outfile << f_tn[i] << " ";
                if (i==nx-1) {
                    outfile << std::endl;
                }
#pragma omp task
                {
                f_tm[i]= f_t[i];
                f_t[i] = f_tn[i];
                }
            }
        }
    }
    else std::cout << "Unable to open file on disk, abort simulation" << std::endl;
    outfile.close();
    
}



