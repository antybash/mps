#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "tensor.h"
#include "DMRG.h"
#include "utilities.h"

int main()
{
    double a = 1.0;
    double b = 1.0; 
    double c = 1.0;

    mpo_heis::initialize(a,b,c);

    std::ofstream outfile;
    std::stringstream in;
    in << "E0_";
    in << std::fixed << std::setprecision(3) << c; 

    double epsilon = 1e-4;
             in << "_1e-4"; 

    in << ".dat";

    std::cout << "Writing to: " << in.str() << std::endl;
    outfile.open(in.str());

    DMRG<18> d16_20(20,epsilon);      d16_20.sweep(2) ;    d16_20.output_lowest_energy(outfile, true); 
                                      d16_20.sweep(2) ;    d16_20.output_lowest_energy(outfile, true); 

    return 0;
}

