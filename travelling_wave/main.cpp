//
//  main.cpp
//  travelling_wave
//
//  Created by Mirco Meazzo on 11/09/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#include <iostream>
#include "omp.h"
#include "support.hpp"

#define NX 100

f wave1(NX);
double ft[NX]= {0};

int main(int argc, const char * argv[]) {
    using namespace std;
    int impulse_pos;
    cout << "Insert the Impulse position.\nSet a value among 1 and " << NX-1 << endl;
    cin >> impulse_pos;
    ft[impulse_pos]=1;
    
//    Preliminar Grid operations
    wave1.f_grid.initialize_grid();
//    wave1.f_grid.show_grid();
    wave1.f_grid.compute_deltax();
    
    
    wave1.init_f(NX,ft);
    wave1.solve_pde(1,15);
    
    
    return 0;
}
