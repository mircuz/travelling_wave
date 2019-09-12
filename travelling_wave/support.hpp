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
    
private:
    double length;
    double* position= new double(nx);
    
protected:
    int nx;
    
public:
    Grid(int n, double len) :
    nx (n), length (len) {
    }
    
    void initialize_grid();
    void show_grid();
    double* compute_deltax();
};



class f : protected Grid{

public:
    double* solve_pde(double c);
    
private:
    double* f_t = new double(nx);
    double* f_tn= new double(nx);   //t+1
    double* f_tm= new double(nx);   //t-1
};


#endif /* support_hpp */
