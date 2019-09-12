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


double *dx;
Grid grid(10,5);


int main(int argc, const char * argv[]) {
    using namespace std;
   
    grid.initialize_grid();
    grid.show_grid();
    dx = grid.compute_deltax();
    
    
    for (int i=0; i<10-1; i++) {
        cout << "dx[" << i << "-" << i+1 << "]=" << *(dx+i) << endl;
    }
    
    
    return 0;
}
