#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "tensor.h"
#include "DMRG.h"
#include "utilities.h"

#ifndef __LENGTH
    #define __LENGTH 20
#endif

#ifndef __MATSIZE
    #define __MATSIZE 20
#endif

#ifndef __EPSILON
    #define __EPSILON 1e-3
#endif

int main()
{
    double a = 1.0;
    double b = 1.0; 
    double c = 1.0;

    mpo_heis::initialize(a,b,c);

    double epsilon = __EPSILON;

    std::stringstream in;
    in << "E0_";
    in << std::fixed << std::setprecision(4) << c; 
    in << "_";
    in << __LENGTH;
    in << "_";
    in << __MATSIZE;
    in << "_";
    in << __EPSILON;

    std::cout << "Writing to: " << in.str() << std::endl;

    DMRG<__LENGTH> dmrg(__MATSIZE, __EPSILON);
    dmrg.set_logs(true, in.str());

    dmrg.sweep(4);
    //dmrg.output_lowest_energy(outfile, true); 
    //dmrg.output_eigenvalue_history(outfile);

    return 0;
}
