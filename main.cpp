#include <iostream>

#include "tensor.h"
#include "DMRG.h"
#include "tests.h"
#include "utilities.h"
#include "linalg.h"

int main()
{
    mpo::initialize();

    tensor<cd,6> psi(vi({1,2,2,2,2,1}));
    psi(vi({0, 0, 0, 0, 0, 0})) = 0;
    psi(vi({0, 1, 0, 0, 0, 0})) = 0;
    psi(vi({0, 0, 1, 0, 0, 0})) = 0;
    psi(vi({0, 1, 1, 0, 0, 0})) = 0;
    psi(vi({0, 0, 0, 1, 0, 0})) = 0; 
    psi(vi({0, 1, 0, 1, 0, 0})) = 1; 
    psi(vi({0, 0, 1, 1, 0, 0})) = 0; 
    psi(vi({0, 1, 1, 1, 0, 0})) = 0;
    psi(vi({0, 0, 0, 0, 1, 0})) = 0;
    psi(vi({0, 1, 0, 0, 1, 0})) = 0;
    psi(vi({0, 0, 1, 0, 1, 0})) = 0;
    psi(vi({0, 1, 1, 0, 1, 0})) = 0;
    psi(vi({0, 0, 0, 1, 1, 0})) = 0; 
    psi(vi({0, 1, 0, 1, 1, 0})) = 0; 
    psi(vi({0, 0, 1, 1, 1, 0})) = 0; 
    psi(vi({0, 1, 1, 1, 1, 0})) = 0;

    auto mps = tensor_to_left_normalized_mps(psi,-1);
    assert(mps.size() == 4);

    std::vector<tensor<cd,4> > mpoHam;
    mpoHam.push_back( mpo::startH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::endH );

    auto eig = DMRG_sweep(mps,mpoHam);
    for(auto i : eig)
        std::cout << i << std::endl;

    output_tensorfull( contract(mps[0],contract(mps[1],contract(mps[2],mps[3],ar::two,ar::zero),ar::two,ar::zero),ar::two,ar::zero) ); 

    return 0;
}

