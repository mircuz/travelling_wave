//
//  support.hpp
//  travelling_wave
//
//  Created by Mirco Meazzo on 11/09/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//


/*------------------ HEADER FILE ------------------*/
#ifndef support_hpp
#define support_hpp

#include <stdio.h>



class Grid {
    
protected:
    int nx;
    int ny;
    
private:
    double length_x, length_y;
    double* positionx= new double(nx);
    double* positiony= new double(ny);
    
public:
    Grid(int nx, int ny, double lenx=1, double leny=1) :
    nx (nx), ny (ny), length_x (lenx), length_y (leny) {
    }
    
    void initialize_grid();
    void show_grid();
    double compute_deltax();
    double compute_deltay();
};



class f {

public:
    int nx;
    int ny;
    Grid f_grid;
    f(int nx_, int ny_) : nx(nx_), ny(ny_), f_grid(nx_,ny_) {   /*Request number of nodes*/
    };
    void solve_pde(double c, double t_limit);
    void init_f( int nx, int ny, double arr[][ny]);
protected:
    double** f_t = new double*[nx];
    double** f_tn= new double*[nx];   //t+1
    double** f_tm= new double*[nx];   //t-1
};


#endif /* support_hpp */
