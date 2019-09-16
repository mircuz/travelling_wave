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

#define NX 10
#define NY 10

f wave1(NX,NY);
double ft[NX][NY]= {0};

int main(int argc, const char * argv[]) {
    using namespace std;
    int impulsex, impulsey;
    cout << "Insert the Impulse position.\nSet x value among 1 and " << NX-2 << " and y value among 1 and " << NY-2 << endl;
    cin >> impulsex >> impulsey;
    ft[impulsex][impulsey]=1;
    
//    Preliminar Grid operations
    wave1.f_grid.initialize_grid();
//    wave1.f_grid.show_grid();
    
    wave1.init_f(NX,NY,ft);
    wave1.solve_pde(1,2);
    
    
    return 0;
}
