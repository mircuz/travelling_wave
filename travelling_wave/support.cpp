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
#include <cmath>



void Grid::initialize_grid(){
    
    double stepx = length_x / nx;
    double stepy = length_y / ny;
#pragma omp parallel for schedule(dynamic)
    for (int i=0; i<nx; i++) {
        positionx[i]=stepx*i;
    }
    for (int j=0; j<ny; j++) {
        positiony[j]=stepy*j;
    }
}



void Grid::show_grid(){
    
    using namespace std;
    cout << "X:" << endl;
    for (int i=0; i<nx; i++) {
        cout << "[" << i << "]= " << this->positionx[i] << endl;
    }
    cout << "Y:" << endl;
    for (int i=0; i<ny; i++) {
        cout << "[" << i << "]= " << this->positiony[i] << endl;
    }
}



double Grid::compute_deltax(){
    
    double deltax;
    deltax = this->positionx[1] - this->positionx[0];
    
    return deltax;
}

double Grid::compute_deltay(){
    double deltay;
    deltay = this->positiony[1] - this->positiony[0];
    
    return deltay;
}




void f::init_f( int nx, int ny, double arr_t[][ny]) {
    
    for (int i=0; i<nx; i++) {
        f_t[i]  = new double[ny];
    }
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            f_t[i][j]= arr_t[i][j];
        }
    }
}



void f::solve_pde(double c, double t_limit) {
    
    for (int i=0; i<nx; i++) {
        f_tn[i] = new double[ny];
        f_tm[i] = new double[ny];
    }
    
    std::ofstream outfile;
    outfile.open("results.txt", std::ios::trunc);
    if (outfile.is_open()) {
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                outfile << f_t[i][j] << "\t";
                if (j==ny-1) {
                    outfile << std::endl;
                }
            }
            if (i==nx-1) {
                outfile << std::endl;
            }
        }
        
        double t=0, dt;
        double** rhs_x = new double*[nx];
        double** rhs_y = new double*[nx];
        for (int i=0; i<nx; i++) {
            rhs_x[i] = new double[ny];
            rhs_y[i] = new double[ny];
        }
        double dx,dy;
        dx = f_grid.compute_deltax();
        dy = f_grid.compute_deltay();
        dt = 1/c * 1/(sqrt((1/(dx*dx) + 1/(dy*dy))));
        std::cout << "dx=" << dx << "\tdy=" << dy << "\tdt=" << dt << std::endl;
        double Cx = c*(dt/dx), Cy = c*(dt/dy);
        std::cout << "Courant number of the simulation set to " << Cx*Cx+Cy*Cy << std::endl;
        
        
        while (t<t_limit) {
//         build rhs
//            Row 0,j
#pragma omp parallel for schedule(dynamic)
            for (int j=0; j<ny; j++) {
                if (j==0) {
                    rhs_x[0][j] = f_t[0][j] + f_t[2][j] -2*f_t[1][j];
                    rhs_y[0][j] = f_t[0][j] + f_t[0][j+2] + -2*f_t[0][j+1];
                }
                if (j==ny-1) {
                    rhs_x[0][j] = f_t[0][j] + f_t[2][j] -2*f_t[1][j];
                    rhs_y[0][j] = f_t[0][j] + f_t[0][j-2] + -2*f_t[0][j-1];
                }
                rhs_x[0][j] = f_t[0][j] + f_t[2][j] -2*f_t[1][j];
                rhs_y[0][j] = f_t[0][j-1] + f_t[0][j+1] -2*f_t[0][j];
            }
//            Row nx-1,j
#pragma omp parallel for schedule(dynamic)
            for (int j=0; j<ny; j++) {
                if (j==0) {
                    rhs_x[nx-1][j] = f_t[nx-1][j] + f_t[nx-3][j] -2*f_t[nx-2][j];
                    rhs_y[nx-1][j] = f_t[nx-1][j] + f_t[nx-1][j+2] + -2*f_t[nx-1][j+1];
                }
                if (j==ny-1) {
                    rhs_x[nx-1][j] = f_t[nx-1][j] + f_t[nx-3][j] -2*f_t[nx-2][j];
                    rhs_y[nx-1][j] = f_t[nx-1][j] + f_t[nx-1][j-2] -2*f_t[nx-1][j-1];
                }
                rhs_x[nx-1][j] = f_t[nx-1][j] + f_t[nx-3][j] -2*f_t[nx-2][j];
                rhs_y[nx-1][j] = f_t[nx-1][j-1] + f_t[nx-1][j+1] -2*f_t[nx-1][j];
            }
//            Col i,0
#pragma omp parallel for schedule(dynamic)
            for (int i=1; i<nx-1; i++) {
                rhs_x[i][0] = f_t[i-1][0] + f_t[i+1][0] - 2*f_t[i][0];
                rhs_y[i][0] = f_t[i][0] + f_t[i][2] - 2*f_t[i][1];
            }
//            Col i,ny-1
#pragma omp parallel for schedule(dynamic)
            for (int i=1; i<nx-1; i++) {
                rhs_x[i][ny-1] = f_t[i-1][ny-1] + f_t[i+1][ny-1] - 2*f_t[i][ny-1];
                rhs_y[i][ny-1] = f_t[i][ny-1] + f_t[i][ny-3] - 2*f_t[i][ny-2];
            }
            
//            Inner points
#pragma omp parallel for schedule(dynamic)
            for (int i=1; i<nx-1; i++) {
                for (int j=1; j<ny-1; j++) {
                    rhs_x[i][j]= f_t[i-1][j] + f_t[i+1][j] - 2*f_t[i][j];
                    rhs_y[i][j]= - 2*f_t[i][j] + f_t[i][j+1] + f_t[i][j-1];
                }
            }
            
            
            if (t==0) {
            //        Solve f at 0+dt
#pragma omp parallel for schedule(dynamic)
                for (int i=0; i<nx; i++) {
                    for (int j=0; j<ny; j++) {
                        f_tn[i][j]= +0.5*Cx*Cx*rhs_x[i][j] +0.5*Cy*Cy*rhs_y[i][j] + f_t[i][j];
                    }
                }
            }
            else {
            //        Solve f at t+dt
#pragma omp parallel for schedule(dynamic)
                for (int i=0; i<nx; i++) {
                    for (int j=0; j<ny; j++) {
                        f_tn[i][j]= Cx*Cx*rhs_x[i][j] + Cy*Cy*rhs_y[i][j] + 2*f_t[i][j]- f_tm[i][j];
                    }
                }
            }
//            wall condition
            for (int i=0; i<nx; i++) {
                f_tn[i][0]=f_tn[i][ny-1]=0;
            }
            for (int j=0; j<ny; j++) {
                f_tn[0][j]=f_tn[nx-1][j]=0;
            }
    
    //        Update f in time
            t=t+dt;
//            Write solution to disk and move it to f_t
            for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                    outfile << f_t[i][j] << "\t";
                    if (j==ny-1) {
                        outfile << std::endl;
                    }
#pragma omp task
                {
                f_tm[i][j]= f_t[i][j];
                f_t[i][j] = f_tn[i][j];
                }
                }
                if (i==nx-1) {
                    outfile << std::endl;
            }
        }
    }
    }
    else std::cout << "Unable to open file on disk, abort simulation" << std::endl;
    outfile.close();
    
}



